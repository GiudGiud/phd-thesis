import glob
import types

import openmoc
import infermc

from discretize import discretize_geometry


statepoints = glob.glob('statepoint.*.h5')

# Initialize a Batchwise for the final statepoint
batchwise = infermc.batchwise.PerfectBatchwise()
batchwise.type = '{}-groups'.format(batchwise.num_fine_groups)
batchwise.sp_start = len(statepoints) - 1
batchwise.cell_mgxslib_filename = 'distribcell'
batchwise.mat_mgxslib_filename = 'material'
batchwise.zcoord = 205.
batchwise.log_level = 'INFO'

# Initialize assembly-wise CMFD mesh
batchwise.cmfd = openmoc.Cmfd()
#batchwise.cmfd.setLatticeStructure(23*17, 23*17)
batchwise.cmfd.setLatticeStructure(23, 23)
batchwise.cmfd.setKNearest(3)
batchwise.with_cmfd = True

# Attach a method to discretize this geometry to the Batchwise instance
batchwise._discretize_geometry = types.MethodType(discretize_geometry, batchwise)

# Assign a clusterizer to the Batchwise
if batchwise.options.clusterizer_type == 'infinite':
    batchwise.clusterizer = infermc.clusterizer.NullClusterizer()
    batchwise.clusterizer._type = 'infinite'
    batchwise.mat_mgxslib_directories.append('../pin-1.6/')
    batchwise.mat_mgxslib_directories.append('../pin-2.4/')
    batchwise.mat_mgxslib_directories.append('../pin-3.1/')
elif batchwise.options.clusterizer_type == 'null':
    batchwise.clusterizer = infermc.clusterizer.NullClusterizer()
elif batchwise.options.clusterizer_type == 'degenerate':
    batchwise.clusterizer = infermc.clusterizer.DegenerateClusterizer()

# Turn off MGXS plotting for speed
batchwise.clusterizer.plot_mgxs = False
batchwise.plot_materials = False
batchwise.plot_cells = False
batchwise.plot_fsrs = False

# Execute OpenMOC simulations over all batches of clustered MGXS libraries
batchwise.execute()
