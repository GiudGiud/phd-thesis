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

        # Extract OpenMC fission rate mesh tally from StatePoint
        curr_fiss = sp.get_tally(name='fission rates')

        # Copy and reshape the NumPy array of mean values
        curr_fiss_mean = copy.deepcopy(curr_fiss.mean)
        curr_fiss_mean.shape = tuple(openmc_mesh.dimension)
        curr_fiss_mean = np.squeeze(curr_fiss_mean)
        curr_fiss_mean = np.fliplr(curr_fiss_mean)

        # Copy and reshape the NumPy array of standard deviation values
        curr_fiss_std_dev = copy.deepcopy(curr_fiss.std_dev)
        curr_fiss_std_dev.shape = tuple(openmc_mesh.dimension)
        curr_fiss_std_dev = np.squeeze(curr_fiss_std_dev)
        curr_fiss_std_dev = np.fliplr(curr_fiss_std_dev)

        # Normalize the capture rates to sum to unity
        curr_fiss_mean /= np.mean(curr_fiss_mean.flat)

        # Set zero capture rates to NaN for transparency in plots
        zero_indices = np.where(curr_fiss_mean < 1E-5)
        curr_fiss_mean[zero_indices] = np.nan
        curr_fiss_std_dev[zero_indices] = np.nan

        # Compute the OpenMC relative error from the non-normalized mean
        curr_rel_err = curr_fiss_std_dev / curr_fiss_mean * 100.

        # Store this batch's relative error
        mean[i, batch-101, :] = curr_fiss_mean.flat
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

plt.title('Max. Fission Rate Error')
plt.grid(True,which="majorminor",ls="-", color='0.65')
plt.xlabel('Batch')
plt.ylabel('Relative Error [%]')
plt.legend(list(directories.values()), loc='upper right')
plt.savefig('fiss-conv-max-assms.png', bbox_inches='tight')
plt.close()

# Create a matplotlib figure for the mean relative error convergence curves
fig = plt.figure()

# Customize and save plot
for i, directory in enumerate(directories):
    plt.loglog(batches, np.nanmean(rel_err[i, :, :], axis=1), linewidth=2)

plt.title('Mean Fission Rate Error')
plt.grid(True,which="majorminor",ls="-", color='0.65')
plt.xlabel('Batch')
plt.ylabel('Relative Error [%]')
plt.legend(list(directories.values()), loc='upper right')
plt.savefig('fiss-conv-mean-assms.png', bbox_inches='tight')
plt.close()


###############################################################################
# MULTIPLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict()
directories['2x2'] = '2x2'
directories['reflector'] = 'reflector'
#directories['full-core'] = 'full core'

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

        # Extract OpenMC fission rate mesh tally from StatePoint
        curr_fiss = sp.get_tally(name='fission rates')

        # Copy and reshape the NumPy array of mean values
        curr_fiss_mean = copy.deepcopy(curr_fiss.mean)
        curr_fiss_mean.shape = tuple(openmc_mesh.dimension)
        curr_fiss_mean = np.squeeze(curr_fiss_mean)
        curr_fiss_mean = np.fliplr(curr_fiss_mean)

        # Copy and reshape the NumPy array of standard deviation values
        curr_fiss_std_dev = copy.deepcopy(curr_fiss.std_dev)
        curr_fiss_std_dev.shape = tuple(openmc_mesh.dimension)
        curr_fiss_std_dev = np.squeeze(curr_fiss_std_dev)
        curr_fiss_std_dev = np.fliplr(curr_fiss_std_dev)

        # Normalize the capture rates to sum to unity
        curr_fiss_mean /= np.mean(curr_fiss_mean.flat)

        # Set zero capture rates to NaN for transparency in plots
        zero_indices = np.where(curr_fiss_mean < 1E-5)
        curr_fiss_mean[zero_indices] = np.nan
        curr_fiss_std_dev[zero_indices] = np.nan

        # Compute the OpenMC relative error from the non-normalized mean
        curr_rel_err = curr_fiss_std_dev / curr_fiss_mean * 100.

        # Store this batch's relative error
        mean[i, batch-101, :] = curr_fiss_mean.flat
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

plt.title('Max. Fission Rate Error')
plt.xlabel('Batch')
plt.ylabel('Relative Error [%]')
plt.legend(list(directories.values()), loc='upper right')
plt.savefig('fiss-conv-max-multi-assms.png', bbox_inches='tight')
plt.close()

# Create a matplotlib figure for the mean relative error convergence curves
fig = plt.figure()

# Customize and save plot
for i, directory in enumerate(directories):
    plt.loglog(batches, np.nanmean(rel_err[i, :, :], axis=1), linewidth=2)

plt.title('Mean Fission Rate Error')
plt.grid(True,which="majorminor",ls="-", color='0.65')
plt.xlabel('Batch')
plt.ylabel('Relative Error [%]')
plt.legend(list(directories.values()), loc='upper right')
plt.savefig('fiss-conv-mean-multi-assms.png', bbox_inches='tight')
plt.close()
