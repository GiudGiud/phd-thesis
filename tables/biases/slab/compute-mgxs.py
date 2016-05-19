import openmc
import openmc.mgxs
from infermc.energy_groups import group_structures

# Load the last statepoint and summary files
sp = openmc.StatePoint('statepoint.100.h5')

# Initialize a fine (70-) group MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(sp.summary.openmc_geometry, by_nuclide=False)
mgxs_lib.energy_groups = group_structures['CASMO']['8-group']
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'cell'
mgxs_lib.correction = None
mgxs_lib.build_library()
mgxs_lib.load_from_statepoint(sp)
mgxs_lib.dump_to_file()
