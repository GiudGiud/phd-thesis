import numpy as np

import openmc
import openmc.opencg_compatible as opencg_compatible
import infermc.beavrs
from infermc.energy_groups import group_structures


####################   User-specified Simulation Parameters  ###################

batches = 10000
inactive = 100
particles = 100000


#########   Exporting to OpenMC geometry.xml and materials.xml Files  ##########

# Write all BEAVRS materials to materials.xml file
infermc.beavrs.make_iso_in_lab()
infermc.beavrs.write_materials_file()

# Extract reflected geometry from InferMC's pre-built assembly Geometries
reflector = infermc.beavrs.build_reflector('Fuel 1.6% enr instr no BAs',
                                           'Fuel 3.1% enr instr 20')
openmc_geometry = opencg_compatible.get_openmc_geometry(reflector)
openmc_geometry.export_to_xml()

##################   Exporting to OpenMC settings.xml File  ###################

# Construct uniform initial source distribution over fissionable zones
lower_left = reflector.bounds[:3]
upper_right = reflector.bounds[3:]

lat_width = (np.array(upper_right) - np.array(lower_left))
lat_width[:2] /= 3.

source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
source.space.only_fissionable = True

settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.statepoint_interval = 20
settings_file.ptables = True
settings_file.output = {'tallies': False}
settings_file.source = source
settings_file.sourcepoint_write = False
settings_file.export_to_xml()


########################   Build OpenMC MGXS Library  #########################

# Get all cells filled with a "fuel" material
mat_cells = openmc_geometry.get_all_material_cells()
fuel_cells = []
for cell in mat_cells:
    if 'fuel' in cell.fill.name.lower():
        fuel_cells.append(cell)

# Initialize a fine (70-) group "distribcell" MGXS Library for OpenMOC
cell_mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
cell_mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
cell_mgxs_lib.mgxs_types = ['total', 'fission', 'capture']
cell_mgxs_lib.domain_type = 'distribcell'
cell_mgxs_lib.domains = fuel_cells
cell_mgxs_lib.correction = None
cell_mgxs_lib.build_library()

# Select the nuclides for the distribcell MGXS
for domain in cell_mgxs_lib.domains:
    for mgxs_type in cell_mgxs_lib.mgxs_types:
        mgxs = cell_mgxs_lib.get_mgxs(domain.id, mgxs_type)
        mgxs.nuclides = [*mgxs.nuclides, 'total']

# Initialize a fine (70-) group "material" MGXS Library for OpenMOC
mat_mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
mat_mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mat_mgxs_lib.mgxs_types = ['total', 'fission', 'capture']
mat_mgxs_lib.domain_type = 'material'
mat_mgxs_lib.correction = None
mat_mgxs_lib.build_library()

# Create a "tallies.xml" file for the mesh tallies
tallies_file = openmc.Tallies()
cell_mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
mat_mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()
