import glob
import types

import sklearn.cluster
import openmc
import infermc
from discretize import discretize_geometry


# Initialize a Batchwise for the final statepoint
batchwise = infermc.batchwise.PerfectBatchwise()
batchwise.type = '{}-groups'.format(batchwise.num_fine_groups)
batchwise.cell_mgxslib_filename = 'distribcell'
batchwise.mat_mgxslib_filename = 'material'
batchwise.zcoord = 195.
batchwise.log_level = 'INFO'
batchwise.with_cmfd = True
batchwise.reference_sp = openmc.StatePoint('../../benchmarks/reflector/statepoint.1000.h5')

# Attach a method to discretize this geometry to the Batchwise instance
batchwise._discretize_geometry = types.MethodType(discretize_geometry, batchwise)

# Initialize a scikit-learn clustering estimator
estimator = sklearn.cluster.AgglomerativeClustering()
estimator.n_clusters = batchwise.options.num_clusters

if batchwise.options.clusterizer_type == 'null':
    batchwise.clusterizer = infermc.clusterizer.NullClusterizer()
elif batchwise.options.clusterizer_type == 'degenerate':
    batchwise.clusterizer = infermc.clusterizer.DegenerateClusterizer()
elif batchwise.options.clusterizer_type == 'combined':
    clusterizer = infermc.clusterizer.CombinedClusterizer()
    batchwise.clusterizer = clusterizer
    batchwise.clusterizer.estimator = estimator
elif batchwise.options.clusterizer_type == 'pinch':
    clusterizer = infermc.clusterizer.PinchClusterizer(6e-6, 'U238', 'capture')
    batchwise.clusterizer = clusterizer
    batchwise.clusterizer.estimator = estimator

# Turn off MGXS plotting for speed
if batchwise.clusterizer:
    batchwise.clusterizer.plot_mgxs = False

batchwise.plot_materials = True
batchwise.plot_cells = False
batchwise.plot_fsrs = True

# Execute OpenMOC simulations over all batches of clustered MGXS libraries
batchwise.execute()
