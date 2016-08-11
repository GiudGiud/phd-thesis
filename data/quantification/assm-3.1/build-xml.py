import numpy as np

import openmc
import openmc.opencg_compatible as opencg_compatible
import infermc.beavrs
from infermc.energy_groups import group_structures


####################   User-specified Simulation Parameters  ###################

batches = 1000
inactive = 100
particles = 1000000


#########   Exporting to OpenMC geometry.xml and materials.xml Files  ##########

# Write all BEAVRS materials to materials.xml file
infermc.beavrs.make_iso_in_lab()
infermc.beavrs.write_materials_file()

# Extract fuel assembly of interest from BEAVRS model
fuel_assembly = infermc.beavrs.find_assembly('Fuel 3.1% enr instr no BAs')
openmc_geometry = opencg_compatible.get_openmc_geometry(fuel_assembly)
openmc_geometry.export_to_xml()


##################   Exporting to OpenMC settings.xml File  ###################

# Construct uniform initial source distribution over fissionable zones
lower_left = fuel_assembly.bounds[:3]
upper_right = fuel_assembly.bounds[3:]
source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
source.space.only_fissionable = True

settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.statepoint_interval = 10
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
cell_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                            'chi', 'absorption', 'capture']
cell_mgxs_lib.domain_type = 'distribcell'
cell_mgxs_lib.domains = fuel_cells
cell_mgxs_lib.correction = None
cell_mgxs_lib.build_library()

# Initialize a fine (70-) group "material" MGXS Library for OpenMOC
mat_mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
mat_mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mat_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                           'chi', 'absorption', 'capture']
mat_mgxs_lib.domain_type = 'material'
mat_mgxs_lib.correction = None
mat_mgxs_lib.build_library()


###################  Create Mesh Tallies for Verification  ####################

# Instantiate a tally Mesh
mesh = openmc.Mesh(name='assembly mesh')
mesh.type = 'regular'
mesh.dimension = [17, 17, 1]
mesh.lower_left = lower_left
mesh.width = (np.array(upper_right) - np.array(lower_left))
mesh.width[:2] /= 17

# Instantiate tally Filter
mesh_filter = openmc.Filter()
mesh_filter.mesh = mesh

# Instantiate energy-integrated fission rate mesh Tally
fission_rates = openmc.Tally(name='fission rates')
fission_rates.filters = [mesh_filter]
fission_rates.scores = ['fission']

# Instantiate energy-wise U-238 capture rate mesh Tally
capture_rates = openmc.Tally(name='u-238 capture')
capture_rates.filters = [mesh_filter]
capture_rates.nuclides = ['U238']
capture_rates.scores = ['absorption', 'fission']

# Create a "tallies.xml" file for the mesh tallies
tallies_file = openmc.Tallies([fission_rates, capture_rates])
cell_mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
mat_mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()
