import numpy as np
import matplotlib 

# Force headless backend for plotting on clusters
matplotlib.use('Agg')

import matplotlib.pylab as pylab
import seaborn as sns

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from infermc.energy_groups import group_structures

# Force non-interactive mode for plotting on clusters
pylab.ioff()

sns.set_style('ticks')


openmoc.log.set_log_level('NORMAL')
opts = openmoc.options.Options()

# Query the user for the number of energy groups
scatter = 'iso-in-lab'
mesh = 32
num_groups = 16

directory = '{}/{}x'.format(scatter, mesh)

# Load the last statepoint and summary files
sp = openmc.StatePoint(directory + '/statepoint.100.h5')

# Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

# Build a coarse group Library from the fine (70-)group Library
coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
condense_lib = mgxs_lib.get_condensed_library(coarse_groups)
condense_lib = condense_lib.get_subdomain_avg_library()

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = get_openmoc_geometry(condense_lib.opencg_geometry)
openmoc.materialize.load_openmc_mgxs_lib(condense_lib, openmoc_geometry)

# Generate tracks
track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()

# Run OpenMOC
solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(opts.num_omp_threads)
solver.setConvergenceThreshold(1E-7)
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()

# Print report of keff and bias with OpenMC
openmoc_keff = solver.getKeff()
openmc_keff = sp.k_combined[0]
bias = (openmoc_keff - openmc_keff) * 1e5

openmoc.log.py_printf('RESULT', 'OpenMOC keff = %f', openmoc_keff)
openmoc.log.py_printf('RESULT', 'OpenMC keff = %f', openmc_keff)
openmoc.log.py_printf('RESULT', 'Bias [pcm] = %1.1f', bias)

# Allocate arrays for FSR-specific data to extract from OpenMOC model
num_fsrs = openmoc_geometry.getNumFSRs()
cell_ids = np.zeros(num_fsrs, dtype=np.int)
centroids = np.zeros(num_fsrs, dtype=np.float)
volumes = np.zeros(num_fsrs, dtype=np.float)
openmoc_fluxes = openmoc.process.get_scalar_fluxes(solver)

# Find the cell IDs, volumes, centroids and fluxes for each FSR
for fsr_id in range(num_fsrs):
    cell = openmoc_geometry.findCellContainingFSR(fsr_id)
    cell_ids[fsr_id] = cell.getId()
    volumes[fsr_id] = solver.getFSRVolume(fsr_id)
    centroids[fsr_id] = cell.getMinX()

# Organize cell IDs, volumes and fluxes in order of increasing centroid
indices = np.argsort(centroids)
centroids = centroids[indices]
cell_ids = cell_ids[indices]
volumes = volumes[indices]
openmoc_fluxes = openmoc_fluxes[indices, ::-1]

# Set centroids to the min/max x values 
centroids[0] = 0.
centroids[-1] = 5.

# Get OpenMC fluxes
tot_fiss_src = 0.
openmc_fluxes = np.zeros((num_fsrs, coarse_groups.num_groups))
for fsr_id, cell_id in enumerate(cell_ids):

    # Get NuFissionXS for cell from MGXS Library
    mgxs = condense_lib.get_mgxs(cell_id, 'nu-fission')

    # Store this cell's flux
    flux_tally = mgxs.tallies['flux']
    openmc_fluxes[fsr_id, :] = flux_tally.get_values().flatten()

    # Increment the total fission source
    nu_fission = mgxs.tallies['nu-fission']
    tot_fiss_src += np.sum(nu_fission.mean)

# Normalize flux to total fission source and volume
openmc_fluxes /= volumes[:,np.newaxis] * tot_fiss_src

# Plot the OpenMOC and OpenMC spatially-varying fluxes
for group in range(coarse_groups.num_groups):
    fig = pylab.figure()
    pylab.plot(centroids, openmoc_fluxes[:,group])
    pylab.plot(centroids, openmc_fluxes[:,group])
    pylab.legend(['OpenMOC', 'OpenMC'], loc='best', fontsize=12)
    pylab.title('Scalar Flux (Group {}/{})'.format(group+1, num_groups), y=1.03, fontsize=16)
    pylab.xlabel('x [cm]', fontsize=12)
    pylab.ylabel('Flux', fontsize=12)
    pylab.grid()
    pylab.savefig('flux-group-{}.png'.format(group+1), bbox_inches='tight')
