import glob
import types

import h5py
import numpy as np

import infermc
from discretize import discretize_geometry


statepoints = glob.glob('statepoint.*.h5')
groups = [2, 8] #, 70]
clusterizer_types = ['infinite', 'null'] #, 'degenerate']

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

print('MAX FISSION RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{1.6\% Assm.}} '
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']['openmoc rel. err.']
        msg += '& {:1.2E} '.format(np.nanmax(np.ravel(bias)))
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
        msg += '& {:1.2E} '.format(np.nanmean(np.ravel(bias)))
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
        msg += '& {:1.2E} '.format(np.nanmax(np.ravel(bias)))
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
        msg += '& {:1.2E} '.format(np.nanmean(np.ravel(bias)))
        f.close()
    msg += '\\\\\n'
msg += '\\\\\n'
print(msg)