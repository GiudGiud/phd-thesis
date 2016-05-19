import numpy as np

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from infermc.energy_groups import group_structures


openmoc.log.set_log_level('RESULT')
opts = openmoc.options.Options()

# Query the user for the number of energy groups
num_groups = int(input('Enter # energy_groups: '))

# Load the last statepoint and summary files
sp = openmc.StatePoint('statepoint.100.h5')

# Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
mgxs_lib = openmc.mgxs.Library.load_from_file()

# Build a coarse group Library from the fine (70-)group Library
coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
condense_lib = mgxs_lib.get_condensed_library(coarse_groups)
condense_lib = condense_lib.get_subdomain_avg_library()

# Create an OpenMOC Geometry from the OpenCG Geometry
openmoc_geometry = get_openmoc_geometry(condense_lib.opencg_geometry)
openmoc.materialize.load_openmc_mgxs_lib(condense_lib, openmoc_geometry)

spacings = np.array([0.1, 0.01, 0.001], dtype=np.float)
angles = np.array([4, 8, 16, 32, 64, 128], dtype=np.int)
keffs = np.zeros((len(angles), len(spacings)), dtype=np.float)

for i, angle in enumerate(angles):
    for j, spacing in enumerate(spacings):
        print('angles = {}, spacing = {}'.format(angle, spacing))

        # Generate tracks
        track_generator = openmoc.TrackGenerator(openmoc_geometry, int(angle), spacing)
        track_generator.setNumThreads(opts.num_omp_threads)
        track_generator.generateTracks()

        # Run OpenMOC
        solver = openmoc.CPUSolver(track_generator)
        solver.setNumThreads(opts.num_omp_threads)
        solver.setConvergenceThreshold(opts.tolerance)
        solver.computeEigenvalue(opts.max_iters)
        keffs[i,j] = solver.getKeff()

# Compute the bias with OpenMC in units of pcm
biases = (keffs - sp.k_combined[0]) * 1e5

print(biases)

# Print table for LaTeX
for i, angle in enumerate(angles):
    print('{} & {:1.1f} & {:1.1f} & {:1.1f} \\\\'.format(angle, biases[i,0], biases[i,1], biases[i,2]))