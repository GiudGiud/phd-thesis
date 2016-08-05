import os
import copy
import openmc
import numpy as np
import matplotlib.pyplot as plt
from _collections import OrderedDict

###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict({'fuel-1.6': '1.6% Enr. (no BPs)',
                           'fuel-3.1':'3.1% Enr. (no BPs)',
                           'fuel-3.1-20BAs': '3.1% Enr. (20 BPs)',
                           '2x2': '2x2 Colorset',
                           'reflector': '2x2 Colorset w/ Reflector',
                           'full-core': 'Full Core'})

for directory in directories:
    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.1000.h5'))
    batches = np.linspace(1, sp.current_batch, sp.current_batch, dtype=np.int)

    # Extract mesh for mesh tallies from OpenMC StatePoint
    openmc_mesh = sp.meshes[10000]

    # Extract OpenMC fission rate mesh tally from StatePoint
    fiss = sp.get_tally(name='fission rates')

    # Copy and reshape the NumPy array of mean values
    fiss_mean = copy.deepcopy(fiss.mean)
    fiss_mean.shape = tuple(openmc_mesh.dimension)
    fiss_mean = np.squeeze(fiss_mean)
    fiss_mean = np.fliplr(fiss_mean)

    # Copy and reshape the NumPy array of standard deviation values
    fiss_std_dev = copy.deepcopy(fiss.std_dev)
    fiss_std_dev.shape = tuple(openmc_mesh.dimension)
    fiss_std_dev = np.squeeze(fiss_std_dev)
    fiss_std_dev = np.fliplr(fiss_std_dev)

    # Normalize the fission rates to sum to unity
    fiss_mean /= np.mean(fiss_mean.flat)

    # Set zero fission rates to NaN for transparency in plots
    zero_indices = np.where(fiss_mean < 1E-3)
    fiss_mean[zero_indices] = np.nan
    fiss_std_dev[zero_indices] = np.nan

    # Compute the OpenMC relative error from the non-normalized mean
    rel_err = fiss_std_dev / fiss_mean * 100.

    # Instantiate a color mapping scheme which makes NaNs transparent
    cmap = plt.get_cmap('jet')
    cmap.set_bad(alpha=0.0)

    # Create a matplotlib figure for mean fission rates
    fig = plt.figure()
    plt.imshow(fiss_mean, interpolation='none', cmap=cmap,
               vmin=np.nanmin(fiss_mean), vmax=np.nanmax(fiss_mean))
    cbar = plt.colorbar(format='%1.3f')
    cbar.ax.ticklabel_format(fontsize=20)
    plt.title('Normalized Fission Rates')
    plt.grid(False)
    plt.axis('off')
    plt.savefig('fiss-mean-' + directory + '.png', bbox_inches='tight')

    # Create a matplotlib figure for rel. err. fission rates
    fig = plt.figure()
    plt.imshow(rel_err, interpolation='none', cmap=cmap,
               vmin=np.nanmin(rel_err), vmax=np.nanmax(rel_err))
    cbar = plt.colorbar(format='%1.2e')
    cbar.ax.ticklabel_format(fontsize=20)
    plt.title('Fission Rate Relative Error [%]')
    plt.grid(False)
    plt.axis('off')
    plt.savefig('fiss-rel-err-' + directory + '.png', bbox_inches='tight')