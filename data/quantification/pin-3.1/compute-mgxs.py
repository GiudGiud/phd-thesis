import glob

import openmc
import openmc.mgxs
from infermc.energy_groups import group_structures

# Load the summary file
statepoints = sorted(glob.glob('statepoint.*.h5'))
sp = openmc.StatePoint(statepoints[-1])

# Initialize a fine (70-) group "material" MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(sp.summary.openmc_geometry, by_nuclide=True)
mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                       'chi', 'absorption', 'capture']
mgxs_lib.domain_type = 'material'
mgxs_lib.correction = None
mgxs_lib.build_library()
mgxs_lib.load_from_statepoint(sp)

# Select the nuclides for the material MGXS
for domain in mgxs_lib.domains:
    for mgxs_type in mgxs_lib.mgxs_types:
        mgxs = mgxs_lib.get_mgxs(domain.id, mgxs_type)
        mgxs.nuclides = [*mgxs.nuclides, 'total']

mgxs_lib.dump_to_file(filename='material')
