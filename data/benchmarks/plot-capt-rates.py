import os
import copy

import h5py
import numpy as np
import matplotlib
from _collections import OrderedDict

import openmc

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict({'assm-1.6': '1.6% Enr. (no BPs)',
                           'assm-3.1':'3.1% Enr. (no BPs)',
                           'assm-3.1-20BPs': '3.1% Enr. (20 BPs)',
                           '2x2': '2x2 Colorset',
                           'reflector': '2x2 Colorset w/ Reflector',
                           'full-core': 'BEAVRS Quarter Core'})

for directory in directories:
    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.1000.h5'))
    batches = np.linspace(1, sp.current_batch, sp.current_batch, dtype=np.int)

    # Extract mesh for mesh tallies from OpenMC StatePoint
    openmc_mesh = sp.meshes[10000]

    # Extract OpenMC capture rate mesh tally from StatePoint
    capt = sp.get_tally(name='u-238 capture')
    absorb = capt.get_slice(scores=['absorption'])
    fission = capt.get_slice(scores=['fission'])
    capt = absorb - fission

    # Normalize the capture rates to sum to unity
    capt_mean = copy.deepcopy(capt.mean)
    capt /= np.mean(np.ravel(capt_mean[capt_mean != 0.]))
    
    # Copy and reshape the NumPy array of mean values
    capt_mean = copy.deepcopy(capt.mean)
    capt_mean.shape = tuple(openmc_mesh.dimension)
    capt_mean = np.squeeze(capt_mean)
    capt_mean = np.fliplr(capt_mean)
    capt_mean = capt_mean[::-1, ::-1]

    # Copy and reshape the NumPy array of standard deviation values
    capt_std_dev = copy.deepcopy(capt.std_dev)
    capt_std_dev.shape = tuple(openmc_mesh.dimension)
    capt_std_dev = np.squeeze(capt_std_dev)
    capt_std_dev = np.fliplr(capt_std_dev)
    capt_std_dev = capt_std_dev[::-1, ::-1]

    # Set zero capture rates to NaN for transparency in plots
    zero_indices = np.where(capt_mean < 1E-5)
    capt_mean[zero_indices] = np.nan
    capt_std_dev[zero_indices] = np.nan

    # Extract the array and make it quarter core symmetric                  
    if 'full-core' in directory:
        full_mean = np.zeros((2*capt_mean.shape[0]-1, 2*capt_mean.shape[0]-1), \
dtype=np.float)
        full_mean[:capt_mean.shape[0], :capt_mean.shape[1]] = capt_mean[:,::-1]
        full_mean[capt_mean.shape[0]-1:, :capt_mean.shape[1]] = capt_mean[::-1,::-1]
        full_mean[:capt_mean.shape[0], capt_mean.shape[1]-1:] = capt_mean[:,:]
        full_mean[capt_mean.shape[0]-1:, capt_mean.shape[1]-1:] = capt_mean[::-1,:]  
        full_mean[capt_mean.shape[0]-1,:] *= 2.
        full_mean[:,capt_mean.shape[1]-1] *= 2.

        full_std_dev = np.zeros((2*capt_std_dev.shape[0]-1, 2*capt_std_dev.shape[0]-1), \
dtype=np.float)
        full_std_dev[:capt_std_dev.shape[0], :capt_std_dev.shape[1]] = capt_std_dev[:,::-1]
        full_std_dev[capt_std_dev.shape[0]-1:, :capt_std_dev.shape[1]] = capt_std_dev[::-1,::-1]
        full_std_dev[:capt_std_dev.shape[0], capt_std_dev.shape[1]-1:] = capt_std_dev[:,:]
        full_std_dev[capt_std_dev.shape[0]-1:, capt_std_dev.shape[1]-1:] = capt_std_dev[::-1,:]  
        full_std_dev[capt_std_dev.shape[0]-1,:] *= 1.26492**2
        full_std_dev[:,capt_std_dev.shape[1]-1] *= 1.26492**2

        capt_mean = full_mean
        capt_std_dev = full_std_dev

    # Compute the OpenMC relative error from the non-normalized mean
    rel_err = capt_std_dev / capt_mean * 100.

    # Instantiate a color mapping scheme which makes NaNs transparent
    cmap = plt.get_cmap('jet')
    cmap.set_bad(alpha=0.0)

    # Create a matplotlib figure for mean capture rates
    fig = plt.figure()
    plt.imshow(capt_mean, interpolation='none', cmap=cmap,
               vmin=np.nanmin(capt_mean), vmax=np.nanmax(capt_mean))
    cbar = plt.colorbar(format='%1.3f')
    cbar.ax.ticklabel_format(fontsize=20)
    plt.title('Normalized U-238 Capture Rates', fontsize=20)
    plt.grid(False)
    plt.axis('off')
    plt.savefig('capt-mean-' + directory.replace('.', '') + '.png', bbox_inches='tight')

    # Create a matplotlib figure for rel. err. fission rates
    fig = plt.figure()
    plt.imshow(rel_err, interpolation='none', cmap=cmap,
               vmin=np.nanmin(rel_err), vmax=np.nanmax(rel_err))
    cbar = plt.colorbar(format='%1.2e')
    cbar.ax.ticklabel_format(fontsize=20)
    plt.title('U-238 Capture Rate Relative Error [%]', fontsize=20)
    plt.grid(False)
    plt.axis('off')
    plt.savefig('capt-rel-err-' + directory.replace('.', '') + '.png', bbox_inches='tight')
