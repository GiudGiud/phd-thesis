"""This script generates Table 5.13 - eigenvalue bias by the number of azimuthal
angles and track spacing, and anisotropic vs. isotropic in lab scattering."""


import numpy as np

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry


openmoc.log.set_log_level('RESULT')
opts = openmoc.options.Options()

spacings = np.array([0.1, 0.01, 0.001], dtype=np.float)
angles = np.array([4, 8, 16, 32, 64, 128, 256, 512], dtype=np.int)
scattering = ['anisotropic', 'iso-in-lab']
keffs = np.zeros((len(scattering), len(angles), len(spacings)), dtype=np.float)
biases = np.zeros((len(scattering), len(angles), len(spacings)), dtype=np.float)

for i, scatter in enumerate(scattering):
    print(scatter)

    for j, angle in enumerate(angles):
        for k, spacing in enumerate(spacings):
            print('angles = {}, spacing = {}'.format(angle, spacing))

            # Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
            directory = '{}/1x/'.format(scatter)
            sp = openmc.StatePoint(directory + 'statepoint.100.h5')
            mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

            # Create an OpenMOC Geometry from the OpenCG Geometry
            openmoc_geometry = get_openmoc_geometry(mgxs_lib.opencg_geometry)
            openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

            # Apply radial/angular discretization mesh
            cells = openmoc_geometry.getAllMaterialCells()
            for cell_id, cell in cells.items():
                cell.setNumSectors(8)
                if cell.getName() in ['fuel', 'water']:
                    cell.setNumRings(5)

            # Generate tracks
            track_generator = openmoc.TrackGenerator(openmoc_geometry, int(angle), spacing)
            track_generator.setNumThreads(opts.num_omp_threads)
            track_generator.generateTracks()

            # Run OpenMOC
            solver = openmoc.CPUSolver(track_generator)
            solver.setNumThreads(opts.num_omp_threads)
            solver.setConvergenceThreshold(1E-7)
            solver.computeEigenvalue(opts.max_iters)
            keffs[i,j,k] = solver.getKeff()

    # Compute the bias with OpenMC in units of pcm for this scattering type
    biases[i, ...] = (keffs[i, ...] - sp.k_combined[0]) * 1e5

print(biases)

# Print table for LaTeX
for j, angle in enumerate(angles):
    row = '{}'.format(angle)
    for i, scatter in enumerate(scattering):
        row += ' &'
        for k, spacing in enumerate(spacings):
            row += ' {:1.0f} &'.format(biases[i,j,k])
    print(row[:-1] + '\\\\')
