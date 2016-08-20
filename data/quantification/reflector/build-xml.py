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

# Select the nuclides for the distribcell MGXS
for domain in cell_mgxs_lib.domains:
    for mgxs_type in cell_mgxs_lib.mgxs_types:
        mgxs = cell_mgxs_lib.get_mgxs(domain.id, mgxs_type)
        mgxs.nuclides = [*mgxs.nuclides, 'total']

# Initialize a fine (70-) group "material" MGXS Library for OpenMOC
mat_mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=False)
mat_mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mat_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                           'chi', 'absorption', 'capture']
mat_mgxs_lib.domain_type = 'material'
mat_mgxs_lib.correction = None
mat_mgxs_lib.build_library()

# Select the nuclides for the material MGXS
for domain in mat_mgxs_lib.domains:
    for mgxs_type in mat_mgxs_lib.mgxs_types:
        mgxs = mat_mgxs_lib.get_mgxs(domain.id, mgxs_type)
        mgxs.nuclides = [*mgxs.nuclides, 'total']


###################  Create Mesh Tallies for Verification  ####################

# Instantiate a tally Mesh
mesh = openmc.Mesh(name='assembly mesh')
mesh.type = 'regular'
mesh.dimension = [34, 34, 1]
mesh.lower_left = [lower_left[0], lower_left[1] + lat_width[1], lower_left[2]]
mesh.width = np.array(lat_width)
mesh.width[:2] /= 17.

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
