import numpy as np

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from infermc.energy_groups import group_structures


#openmoc.log.set_log_level('RESULT')
opts = openmoc.options.Options()

groups = [1, 2, 4, 8, 16, 25, 40, 70]
scatter = 'iso-in-lab'
mesh = [1, 2, 4, 8, 16]
#mesh = [1, 2, 4, 8, 16, 32, 64]
keffs = np.zeros((2, len(groups), len(mesh)), dtype=np.float)
biases = np.zeros((2, len(groups), len(mesh)), dtype=np.float)


###############################################################################
#                 Eigenvalue Calculation w/o SPH Factors
###############################################################################

print('Without SPH')
for j, num_groups in enumerate(groups):
    for k, num_mesh in enumerate(mesh):
        print('# groups = {}, # mesh = {}'.format(num_groups, num_mesh))

        # Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
        directory = '{}/{}x/'.format(scatter, num_mesh)
        sp = openmc.StatePoint(directory + 'statepoint.100.h5')
        mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

        # Build a coarse group Library from the fine (70-)group Library
        coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
        condense_lib = mgxs_lib.get_condensed_library(coarse_groups)

        # Create an OpenMOC Geometry from the OpenCG Geometry
        openmoc_geometry = get_openmoc_geometry(condense_lib.opencg_geometry)
        openmoc.materialize.load_openmc_mgxs_lib(condense_lib, openmoc_geometry)

        # Generate tracks
        track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
        track_generator.setNumThreads(opts.num_omp_threads)
        track_generator.generateTracks()

        # Instantiate a Solver
        solver = openmoc.CPUSolver(track_generator)
        solver.setNumThreads(opts.num_omp_threads)
        solver.setConvergenceThreshold(1E-7)

        # Run OpenMOC
        solver.computeEigenvalue(opts.max_iters)
        keffs[0,j,k] = solver.getKeff()


###############################################################################
#                 Eigenvalue Calculation w/ SPH Factors
###############################################################################

print('With SPH')
for j, num_groups in enumerate(groups):
    for k, num_mesh in enumerate(mesh):
        print('# groups = {}, # mesh = {}'.format(num_groups, num_mesh))

        # Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
        directory = '{}/{}x/'.format(scatter, num_mesh)
        sp = openmc.StatePoint(directory + 'statepoint.100.h5')
        mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

        # Build a coarse group Library from the fine (70-)group Library
        coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
        condense_lib = mgxs_lib.get_condensed_library(coarse_groups)

        # Create an OpenMOC Geometry from the OpenCG Geometry
        openmoc_geometry = get_openmoc_geometry(condense_lib.opencg_geometry)
        openmoc.materialize.load_openmc_mgxs_lib(condense_lib, openmoc_geometry)

        # Compute SPH factors
        sph, sph_mgxs_lib, sph_indices = \
            openmoc.materialize.compute_sph_factors(
                condense_lib, num_azim=128, azim_spacing=0.05, sph_tol=1E-7,
                num_threads=opts.num_omp_threads, max_sph_iters=50)

        # Load the SPH-corrected MGXS library data
        materials = \
            openmoc.materialize.load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

        # Generate tracks
        track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.01)
        track_generator.setNumThreads(opts.num_omp_threads)
        track_generator.generateTracks()

        # Instantiate a Solver
        solver = openmoc.CPUSolver(track_generator)
        solver.setNumThreads(opts.num_omp_threads)
        solver.setConvergenceThreshold(1E-7)

        # Run OpenMOC
        solver.computeEigenvalue(opts.max_iters)
        keffs[1,j,k] = solver.getKeff()

# Compute the bias with OpenMC in units of pcm for this scattering type
biases = (keffs - sp.k_combined[0]) * 1e5

print(biases)

# Print w/o SPH table for LaTeX
print('Without SPH')
for i, num_groups in enumerate(groups):
    row = '{} &'.format(num_groups)
    for j, num_mesh in enumerate(mesh):
        row += ' {:1.0f} &'.format(biases[0,i,j])
    print(row[:-1] + '\\\\')

# Print w/ SPH table for LaTeX
print('With SPH')
for i, num_groups in enumerate(groups):
    row = '{} &'.format(num_groups)
    for j, num_mesh in enumerate(mesh):
        row += ' {:1.0f} &'.format(biases[1,i,j])
    print(row[:-1] + '\\\\')
