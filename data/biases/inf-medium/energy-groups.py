import numpy as np

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from infermc.energy_groups import group_structures


openmoc.log.set_log_level('RESULT')
opts = openmoc.options.Options()

groups = [1, 2, 4, 8, 16, 25, 40, 70]

directory = str(input('directory: '))

# Load the last statepoint and summary files
sp = openmc.StatePoint(directory + '/' + 'statepoint.100.h5')

# Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = get_openmoc_geometry(mgxs_lib.opencg_geometry)
openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

# Generate tracks
track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()

# Instantiate a Solver
solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(opts.num_omp_threads)
solver.setConvergenceThreshold(1E-7)

keffs = np.zeros((len(groups),), dtype=np.float)

for i, num_groups in enumerate(groups):
    print('groups = {}'.format(num_groups))

    # Build a coarse group Library from the fine (70-)group Library
    coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
    condense_lib = mgxs_lib.get_condensed_library(coarse_groups)
    openmoc.materialize.load_openmc_mgxs_lib(condense_lib, openmoc_geometry)

    # Run OpenMOC
    solver.computeEigenvalue(opts.max_iters)
    keffs[i] = solver.getKeff()

# Compute the bias with OpenMC in units of pcm
biases = (keffs - sp.k_combined[0]) * 1e5

print(biases)

# Print table for LaTeX
for i, num_groups in enumerate(groups):
    print('{} & {:1.1f} \\\\'.format(num_groups, biases[i]))