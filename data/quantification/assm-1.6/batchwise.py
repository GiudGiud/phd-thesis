import glob
import types
import infermc
from discretize import discretize_geometry

statepoints = glob.glob('statepoint.*.h5')

# FIXME: Loop over clusterizers then groups
groups = [2, 8, 70]
clusterizer_types = ['infinite', 'null', 'degenerate']

for clusterizer_type in clusterizer_types:
    # Initialize a Batchwise to iterate over all statepoints
    batchwise = infermc.batchwise.PerfectBatchwise()
    batchwise.sp_start = len(statepoints) - 1
    batchwise.cell_mgxslib_filename = 'distribcell'
    batchwise.mat_mgxslib_filename = 'material'
    batchwise.zcoord = 205.
    batchwise.log_level = 'INFO'

    # FIXME: Add report to add to table

    # Attach a method to discretize this geometry to the Batchwise instance
    batchwise._discretize_geometry = types.MethodType(discretize_geometry, batchwise)

    if clusterizer_type == 'infinite':
        clusterizer = infermc.clusterizer.NullClusterizer()
        batchwise.clusterizer = clusterizer
        batchwise.clusterizer._type = 'infinite'
        batchwise.mat_mgxslib_directories.append('../pin-1.6/')
    elif clusterizer_type == 'null':
        batchwise.clusterizer = infermc.clusterizer.NullClusterizer()
    elif clusterizer_type == 'degenerate':
        batchwise.clusterizer = infermc.clusterizer.DegenerateClusterizer()

    # Turn off MGXS plotting for speed
    batchwise.clusterizer.plot_mgxs = False

    # Execute OpenMOC simulations over all batches of clustered MGXS libraries
    for num_groups in groups:
        batchwise.num_fine_groups = num_groups
        batchwise.type = '{}-groups'.format(num_groups)
        batchwise.execute()
        batchwise._sp_counter -= 1
