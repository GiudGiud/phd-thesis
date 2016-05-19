import numpy as np

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from infermc.energy_groups import group_structures


def subdivide(cell, num_mesh):
    """Subdivide a cell with a uniform mesh along the x-axis."""

    # Get the Cell's min/max boundaries
    min_x = cell.getMinX()
    max_x = cell.getMaxX()

    # Compute the spacing for the new mesh in the Cell
    delta_x = (max_x - min_x) / num_mesh
    clones = []

    # Subdivide the original cell
    for mesh in range(num_mesh):
        clone = cell.clone()
        clones.append(clone)

        # Adjust X-Planes in the Cell
        surfaces = clone.getSurfaces()
        for surface_id in surfaces:
            surface = surfaces[surface_id]._surface
            if surface.getSurfaceType() == openmoc.XPLANE:
                surface = openmoc.castSurfaceToXPlane(surface)
                if surfaces[surface_id]._halfspace == +1:
                    surface.setX(min_x + i * delta_x)
                else:
                    surface.setX(min_x + (i + 1) * delta_x)

    return clones


openmoc.log.set_log_level('NORMAL')
opts = openmoc.options.Options()

num_mesh = [1, 2, 4, 8, 16, 32, 64]

scattering = str(input('scattering: '))
directory = '{}/1x'.format(scattering)

# Load the last statepoint and summary files
sp = openmc.StatePoint(directory + '/' + 'statepoint.100.h5')

# Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

keffs = np.zeros((len(num_mesh),), dtype=np.float)

for i, mesh in enumerate(num_mesh):

    # Create an OpenMOC Geometry from the OpenCG Geometry
    openmoc_geometry = get_openmoc_geometry(mgxs_lib.opencg_geometry)
    openmoc.materialize.load_openmc_mgxs_lib(mgxs_lib, openmoc_geometry)

    # Discretize the geometry
    '''
    root = openmoc_geometry.getRootUniverse()
    cells = openmoc_geometry.getAllMaterialCells()
    for cell_id in cells:
        print('subdividing cell: {}'.format(cell_id))
        clones = subdivide(cells[cell_id], mesh)
        root.removeCell(cells[cell_id])
        for clone in clones:
            root.addCell(clone)
    '''

    # Generate tracks
    track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
    track_generator.setNumThreads(opts.num_omp_threads)
    track_generator.generateTracks()

    # Instantiate a Solver
    solver = openmoc.CPUSolver(track_generator)
    solver.setNumThreads(opts.num_omp_threads)
    solver.setConvergenceThreshold(1E-7)
    solver.printTimerReport()

    # Run OpenMOC
    solver.computeEigenvalue(opts.max_iters)
    keffs[i] = solver.getKeff()

# Compute the bias with OpenMC in units of pcm
biases = (keffs - sp.k_combined[0]) * 1e5

print(biases)

# Print table for LaTeX
for i, mesh in enumerate(num_mesh):
    print('{} & {:1.1f} \\\\'.format(mesh, biases[i]))