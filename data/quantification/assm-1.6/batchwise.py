import types
import infermc
from discretize import discretize_geometry

# Initialize a Batchwise to iterate over all statepoints
# This requires the --batchwise={'perfect', 'adaptive'} command line argument
batchwise = infermc.batchwise.Batchwise.get_batchwise()
batchwise.cell_mgxslib_filename = 'distribcell'
batchwise.mat_mgxslib_filename = 'material'
batchwise.zcoord = 205.
batchwise.log_level = 'INFO'

# Attach a method to discretize this geometry to the Batchwise instance
batchwise._discretize_geometry = types.MethodType(discretize_geometry, batchwise)

# If a basic clusterizer type ('null', 'degenerate') was not specified
# as a command line argument initialize a clusterizer with scikit-learn
if batchwise.options.clusterizer_type == 'infinite':
    clusterizer = infermc.clusterizer.NullClusterizer()
    batchwise.clusterizer = clusterizer
    batchwise.clusterizer._type = 'infinite'
    batchwise.mat_mgxslib_directories.append('../pin-1.6/')

# Turn off MGXS plotting for speed
if batchwise.clusterizer:
    batchwise.clusterizer.plot_mgxs = False

# Execute OpenMOC simulations over all batches of clustered MGXS libraries
batchwise.execute()
