import os
import copy
import openmc
import numpy as np
from _collections import OrderedDict

import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict()
directories['fuel-1.6'] = '1.6% Enr. (no BPs)'
directories['fuel-3.1'] = '3.1% Enr. (no BPs)'
directories['fuel-3.1-20BAs'] = '3.1% Enr. (20 BPs)'

batches = np.linspace(101, 1000, 1001-101, dtype=np.int)
mean = np.zeros((3, len(batches), 17*17), dtype=np.float)
rel_err = np.zeros((3, len(batches), 17*17), dtype=np.float)

for i, directory in enumerate(directories):
    print(directory)
    for batch in batches:
        print(batch)
        sp = openmc.StatePoint(os.path.join(directory, 'statepoint.{0:04}.h5'.format(batch)))

        # Extract mesh for mesh tallies from OpenMC StatePoint
        openmc_mesh = sp.meshes[10000]

        # Extract OpenMC capture rate mesh tally from StatePoint
        curr_capt = sp.get_tally(name='u-238 capture')
        absorb = curr_capt.get_slice(scores=['absorption'])
        fission = curr_capt.get_slice(scores=['fission'])
        curr_capt = absorb - fission

        # Copy and reshape the NumPy array of mean values
        curr_capt_mean = copy.deepcopy(curr_capt.mean)
        curr_capt_mean.shape = tuple(openmc_mesh.dimension)
        curr_capt_mean = np.squeeze(curr_capt_mean)
        curr_capt_mean = np.fliplr(curr_capt_mean)

        # Copy and reshape the NumPy array of standard deviation values
        curr_capt_std_dev = copy.deepcopy(curr_capt.std_dev)
        curr_capt_std_dev.shape = tuple(openmc_mesh.dimension)
        curr_capt_std_dev = np.squeeze(curr_capt_std_dev)
        curr_capt_std_dev = np.fliplr(curr_capt_std_dev)

        # Normalize the capture rates to sum to unity
        curr_capt_mean /= np.mean(curr_capt_mean.flat)

        # Set zero capture rates to NaN for transparency in plots
        zero_indices = np.where(curr_capt_mean < 1E-5)
        curr_capt_mean[zero_indices] = np.nan
        curr_capt_std_dev[zero_indices] = np.nan

        # Compute the OpenMC relative error from the non-normalized mean
        curr_rel_err = curr_capt_std_dev / curr_capt_mean * 100.

        # Store this batch's relative error
        mean[i, batch-101, :] = curr_capt_mean.flat
        rel_err[i, batch-101, :] = curr_rel_err.flat

# Set the first active batch's uncertainties to NaN
rel_err[:, 0, :] = np.nan

# Create a matplotlib figure for the max relative error convergence curves
fig = plt.figure()

# Customize and save plot
for i, directory in enumerate(directories):
    who = np.nanmax(rel_err[i, :, :], axis=1)
    print(who)
    plt.loglog(batches, np.nanmax(rel_err[i, :, :], axis=1), linewidth=2)

plt.title('Max. U-238 Capture Rate Error')
plt.grid()
plt.xlabel('Batch')
plt.ylabel('Relative Error [%]')
plt.legend(list(directories.values()), loc='upper right')
plt.savefig('capt-conv-max-assms.png', bbox_inches='tight')
plt.close()

# Create a matplotlib figure for the mean relative error convergence curves
fig = plt.figure()

# Customize and save plot
for i, directory in enumerate(directories):
    plt.loglog(batches, np.nanmean(rel_err[i, :, :], axis=1), linewidth=2)

plt.title('Mean U-238 Capture Rate Error')
plt.grid()
plt.xlabel('Batch')
plt.ylabel('Relative Error [%]')
plt.legend(list(directories.values()), loc='upper right')
plt.savefig('capt-conv-mean-assms.png', bbox_inches='tight')
plt.close()
