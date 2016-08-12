import glob

import openmc
import openmc.mgxs
from infermc.energy_groups import group_structures

# Load the summary file
sp = openmc.StatePoint(glob.glob('statepoint.*.h5')[-1])

# Initialize a fine (70-) group "material" MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(sp.summary.openmc_geometry, by_nuclide=True)
mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                       'chi', 'absorption', 'capture']
mgxs_lib.domain_type = 'material'
mgxs_lib.correction = None
mgxs_lib.build_library()
mgxs_lib.load_from_statepoint(sp)
mgxs_lib.dump_to_file(filename='material')