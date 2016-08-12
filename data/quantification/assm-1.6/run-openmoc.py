import glob
import types
import os

import h5py
import numpy as np
import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import infermc
from discretize import discretize_geometry


statepoints = glob.glob('statepoint.*.h5')
groups = [2, 8, 70]
clusterizer_types = ['infinite', 'null', 'degenerate']

# Loop over 'infinite', 'null' and 'degenerate' clusterizers
for clusterizer_type in clusterizer_types:

    # Initialize a Batchwise for the final statepoint
    batchwise = infermc.batchwise.PerfectBatchwise()
    batchwise.sp_start = len(statepoints) - 1
    batchwise.cell_mgxslib_filename = 'distribcell'
    batchwise.mat_mgxslib_filename = 'material'
    batchwise.zcoord = 205.
    batchwise.log_level = 'INFO'

    # Attach a method to discretize this geometry to the Batchwise instance
    batchwise._discretize_geometry = types.MethodType(discretize_geometry, batchwise)

    # Assign a clusterizer to the Batchwise
    if clusterizer_type == 'infinite':
        batchwise.clusterizer = infermc.clusterizer.NullClusterizer()
        batchwise.clusterizer._type = 'infinite'
        batchwise.mat_mgxslib_directories.append('../pin-1.6/')
    elif clusterizer_type == 'null':
        batchwise.clusterizer = infermc.clusterizer.NullClusterizer()
    elif clusterizer_type == 'degenerate':
        batchwise.clusterizer = infermc.clusterizer.DegenerateClusterizer()

    # Turn off MGXS plotting for speed
    batchwise.clusterizer.plot_mgxs = False
    batchwise.plot_materials = False

    # Execute OpenMOC simulations over all batches of clustered MGXS libraries
    for num_groups in groups:
        batchwise.num_fine_groups = num_groups
        batchwise.type = '{}-groups'.format(num_groups)
        batchwise.execute()
        batchwise._sp_counter -= 1

print('EIGENVALUE BIAS')
msg = '\multirow{3}{*}{\parbox{2.5cm}{1.6\% Assm.}} '
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['keff']['bias']
        msg += '& {:.0f} '.format(bias[0])
        f.close()
    msg += '\\\\\n'
msg += '\\\\\n'
print(msg)

# Track the minimum and maximum fission and capture rate errors
min_fiss = +np.inf
max_fiss = -np.inf
min_capt = +np.inf
max_capt = -np.inf

print('MAX FISSION RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{1.6\% Assm.}} '
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']['openmoc rel. err.']
        bias = np.nanmax(np.ravel(bias))
        msg += '& {:1.2E} '.format(bias)
        max_fiss = max(max_fiss, bias)
        min_fiss = min(min_fiss, bias)
        f.close()
    msg += '\\\\\n'
msg += '\\\\\n'
print(msg)

print('MEAN FISSION RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{1.6\% Assm.}} '
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']['openmoc rel. err.']
        bias = np.nanmean(np.ravel(bias))
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
msg += '\\\\\n'
print(msg)

print('MAX CAPTURE RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{1.6\% Assm.}} '
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture']['openmoc rel. err.']
        bias = np.nanmax(np.ravel(bias))
        msg += '& {:1.2E} '.format(bias)
        max_capt = max(max_capt, bias)
        min_capt = min(min_capt, bias)
        f.close()
    msg += '\\\\\n'
msg += '\\\\\n'
print(msg)

print('MEAN CAPTURE RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{1.6\% Assm.}} '
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture']['openmoc rel. err.']
        bias = np.nanmean(np.ravel(bias))
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
msg += '\\\\\n'
print(msg)

if not os.path.exists('plots'):
    os.makedirs('plots')

# Plot fission rate error heat maps
for clusterizer_type in clusterizer_types:
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']

        # Plot a colormap of the fission rate percent rel. err.
        fig = plt.figure()
        cmap = plt.get_cmap('jet')
        cmap.set_bad(alpha=0.0)
        plt.imshow(
            bias['openmoc rel. err.'][-1, ...], interpolation='none', cmap=cmap)
#            vmin=min_fiss, vmax=max_fiss, cmap=cmap)
        plt.title('OpenMOC Rel. Err. [%]', fontsize=12)
        cbar = plt.colorbar()
        cbar.ax.ticklabel_format(fontsize=20)
        plt.grid(False)
        plt.axis('off')
        filename = 'plots/{}-fiss-err-{}.png'.format(clusterizer_type, num_groups)
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)

        f.close()

# Plot U-238 capture rate error heat maps
for clusterizer_type in clusterizer_types:
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture']

        # Plot a colormap of the fission rate percent rel. err.
        fig = plt.figure()
        cmap = plt.get_cmap('jet')
        cmap.set_bad(alpha=0.0)
        plt.imshow(
            bias['openmoc rel. err.'][-1, ...], interpolation='none', cmap=cmap)
#            vmin=min_capt, vmax=max_capt, cmap=cmap)
        plt.title('OpenMOC Rel. Err. [%]', fontsize=12)
        cbar = plt.colorbar()
        cbar.ax.ticklabel_format(fontsize=20)
        plt.grid(False)
        plt.axis('off')
        filename = 'plots/{}-capt-err-{}.png'.format(clusterizer_type, num_groups)
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)

        f.close()
