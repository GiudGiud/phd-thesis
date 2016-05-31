import openmc
import openmc.mgxs
from infermc.energy_groups import group_structures

scattering = ['anisotropic', 'transport', 'iso-in-lab']
mesh = [1, 2, 4, 8, 16]
#mesh = [1, 2, 4, 8, 16, 32]

for scatter in scattering:
    print(scatter)
    for num_mesh in mesh:
        print('# mesh: {}'.format(num_mesh))

        directory = '{}/{}x/'.format(scatter, num_mesh)

        # Load the last statepoint and summary files
        sp = openmc.StatePoint(directory + 'statepoint.100.h5')

        # Initialize a fine (70-) group MGXS Library for OpenMOC
        mgxs_lib = openmc.mgxs.Library(sp.summary.openmc_geometry, by_nuclide=False)
        mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
        if scatter == 'transport':
            mgxs_lib.mgxs_types = ['nu-transport', 'nu-fission', 'nu-scatter matrix', 'chi']
            mgxs_lib.correction = 'P0'
        else:
            mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
            mgxs_lib.correction = None
        mgxs_lib.domain_type = 'cell'
        mgxs_lib.build_library()
        mgxs_lib.load_from_statepoint(sp)
        mgxs_lib.dump_to_file(directory=directory)
