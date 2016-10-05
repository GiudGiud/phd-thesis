import glob
import types

import sklearn.cluster
import sklearn.ensemble
import sklearn.mixture

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
batchwise.reference_sp = openmc.StatePoint('../../benchmarks/assm-3.1/statepoint.1000.h5')

# Attach a method to discretize this geometry to the Batchwise instance
batchwise._discretize_geometry = types.MethodType(discretize_geometry, batchwise)

# Initialize a scikit-learn clustering estimator
if batchwise.options.predictor == 'agglomerative':
    estimator = sklearn.cluster.AgglomerativeClustering()
    estimator.n_clusters = batchwise.options.num_clusters
elif batchwise.options.predictor == 'kmeans':  
    estimator = sklearn.cluster.KMeans()
    estimator.n_clusters = batchwise.options.num_clusters
elif batchwise.options.predictor == 'dbscan':
    estimator = sklearn.cluster.DBSCAN(min_samples=4, eps=0.3, leaf_size=1)
elif batchwise.options.predictor == 'birch':
    estimator = sklearn.cluster.Birch(
        n_clusters=batchwise.options.num_clusters)
elif batchwise.options.predictor == 'gmm':
    estimator = sklearn.mixture.GMM(
        n_components=batchwise.options.num_clusters)
elif batchwise.options.predictor == 'dpgmm':
    estimator = sklearn.mixture.DPGMM(
        n_components=batchwise.options.num_clusters)

# Initialize a Clusterizer
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

#batchwise.clusterizer.litmus = infermc.litmus.DegenerateLitmus()

# FEATURE TRANSFORMATION
if batchwise.options.transformer == 'pca':
    transformer = infermc.transformer.DecomposedTransformer()
    transformer.decomposer = \
        sklearn.decomposition.KernelPCA(n_components=2, kernel='rbf')
    clusterizer.transformer = transformer
# FEATURE REGRESSION
elif batchwise.options.transformer == 'forest':
    regressor = infermc.transformer.TreeBasedRegressor()
    regressor.regressor = sklearn.ensemble.RandomForestRegressor(
        n_estimators=100, max_depth=3, min_samples_leaf=20)
    clusterizer.transformer = regressor  

# Turn off MGXS plotting for speed
if batchwise.clusterizer:
    batchwise.clusterizer.plot_mgxs = False

batchwise.plot_materials = True
batchwise.plot_cells = False
batchwise.plot_fsrs = True

# Execute OpenMOC simulations over all batches of clustered MGXS libraries
batchwise.execute()
