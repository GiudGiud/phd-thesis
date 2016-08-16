import openmc
import openmc.mgxs
from infermc.energy_groups import group_structures

# Load the summary file
su = openmc.Summary('summary.h5')

# Get all cells filled with a "fuel" material
mat_cells = su.openmc_geometry.get_all_material_cells()
fuel_cells = []
for cell in mat_cells:
    if 'fuel' in cell.fill.name.lower():
        fuel_cells.append(cell)

# Initialize a fine (40-) group "distribcell" MGXS Library for OpenMOC
cell_mgxs_lib = openmc.mgxs.Library(su.openmc_geometry, by_nuclide=False)
cell_mgxs_lib.energy_groups = group_structures['CASMO']['25-group']
cell_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                            'chi', 'absorption', 'capture']
cell_mgxs_lib.domain_type = 'distribcell'
cell_mgxs_lib.domains = fuel_cells
cell_mgxs_lib.correction = None
cell_mgxs_lib.build_library()
cell_mgxs_lib.dump_to_file(filename='distribcell')

# Initialize a fine (40-) group "material" MGXS Library for OpenMOC
mat_mgxs_lib = openmc.mgxs.Library(su.openmc_geometry, by_nuclide=False)
mat_mgxs_lib.energy_groups = group_structures['CASMO']['25-group']
mat_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                           'chi', 'absorption', 'capture']
mat_mgxs_lib.domain_type = 'material'
mat_mgxs_lib.correction = None
mat_mgxs_lib.build_library()
mat_mgxs_lib.dump_to_file(filename='material')
