import openmoc
import openmoc.plotter
import openmoc.opencg_compatible
import openmc.mgxs
from infermc.energy_groups import group_structures


opts = openmoc.options.Options()
openmoc.log.set_log_level('NORMAL')

# Load the summary file
su = openmc.Summary('summary.h5')

# Initialize a fine (40-) group "material" MGXS Library for OpenMOC
mat_mgxs_lib = openmc.mgxs.Library(su.openmc_geometry, by_nuclide=True)
mat_mgxs_lib.energy_groups = group_structures['CASMO']['40-group']
mat_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                           'chi', 'absorption', 'capture']
mat_mgxs_lib.domain_type = 'material'
mat_mgxs_lib.correction = None
mat_mgxs_lib.build_library()

# Initialize 40-group OpenMC "material" and MGXS library
sp = openmc.StatePoint('statepoint.0800.h5')
mat_mgxs_lib.load_from_statepoint(sp)

# Condense MGXS library to 2-group structure
coarse_groups = openmc.mgxs.EnergyGroups()
coarse_groups.group_edges = [0., 0.625e-6, 20.]
mat_mgxs_lib = mat_mgxs_lib.get_condensed_library(coarse_groups)

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = \
    openmoc.opencg_compatible.get_openmoc_geometry(mat_mgxs_lib.opencg_geometry)

# Load cross section data
openmoc_materials = \
    openmoc.materialize.load_openmc_mgxs_lib(mat_mgxs_lib, openmoc_geometry)

# Discretize geometry
#from discretize_beavrs import discretize_beavrs
#discretize_beavrs(mat_mgxs_lib, openmoc_geometry)

# Plot the geometry
#openmoc.plotter.plot_cells(openmoc_geometry, zcoord=205.,
#                           gridsize=1000, library='pil')
#openmoc.plotter.plot_materials(openmoc_geometry, zcoord=205.,
#                               gridsize=1000, library='pil')

# Initialize pin-wise CMFD mesh
cmfd = openmoc.Cmfd()
cmfd.setLatticeStructure(23, 23)
cmfd.setKNearest(3)
openmoc_geometry.setCmfd(cmfd)

# Initialize an OpenMOC TrackGenerator and Solver
track_generator = openmoc.TrackGenerator(
    openmoc_geometry, opts.num_azim, opts.azim_spacing)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.setZCoord(205.)
track_generator.generateTracks()

#openmoc.plotter.plot_cmfd_cells(openmoc_geometry, cmfd, gridsize=500)
#openmoc.plotter.plot_flat_source_regions(openmoc_geometry, gridsize=500)

# Initialize an OpenMOC Solver
solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(opts.tolerance)
solver.setNumThreads(opts.num_omp_threads)

# Run an eigenvalue calulation with the MGXS from OpenMC
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()

#openmoc.plotter.plot_residuals(openmoc_geometry)

# Report the bias on the eigenvalue
keff_bias = (solver.getKeff() - sp.k_combined[0]) * 1e5
print('keff bias [pcm]: {:.1f}'.format(keff_bias))
