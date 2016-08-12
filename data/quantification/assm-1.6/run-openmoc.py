import glob
import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from discretize import discretize_geometry
from infermc.energy_groups import group_structures
import infermc


openmoc.log.set_log_level('NORMAL')
opts = openmoc.options.Options()
num_groups = 2

# Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
sp = openmc.StatePoint(glob.glob('statepoint.*.h5')[-1])
cell_mgxs_lib = openmc.mgxs.Library.load_from_file(filename='distribcell')
mat_mgxs_lib = openmc.mgxs.Library.load_from_file(filename='material')
cell_mgxs_lib.load_from_statepoint(sp)
mat_mgxs_lib.load_from_statepoint(sp)

#cell_mgxs_lib = cell_mgxs_lib.get_subdomain_avg_library()

# Build a coarse group Library from the fine (70-)group Library
coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
cell_mgxs_lib = cell_mgxs_lib.get_condensed_library(coarse_groups)
mat_mgxs_lib = mat_mgxs_lib.get_condensed_library(coarse_groups)


clusterizer = infermc.clusterizer.NullClusterizer()

# Assign the MGXS library to the Clusterizer
clusterizer.cell_mgxslib = cell_mgxs_lib


# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = get_openmoc_geometry(cell_mgxs_lib.opencg_geometry)
openmoc.materialize.load_openmc_mgxs_lib(mat_mgxs_lib, openmoc_geometry)
openmoc.materialize.load_openmc_mgxs_lib(cell_mgxs_lib, openmoc_geometry)

# FIXME: Rev your engines for a little discretization....
discretize_geometry(mat_mgxs_lib, openmoc_geometry)

# Initialize CMFD
cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(17,17)
cmfd.setKNearest(3)
openmoc_geometry.setCmfd(cmfd)

# Generate tracks
track_generator = openmoc.TrackGenerator(openmoc_geometry, 32, 0.01)
track_generator.setZCoord(205.0)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()

# FIXME: Attach clusterizer - infinite, null and degenerate

# FIXME: Store keff bias
# FIXME: Compute and store fission rate errors
#  FIXME: Compute and store U-238 capture rate errors
