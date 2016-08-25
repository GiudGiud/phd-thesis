import numpy as np
import openmc
import openmc.opencg_compatible as opencg_compatible
import infermc.beavrs
from infermc.energy_groups import group_structures


####################   User-specified Simulation Parameters  ###################

batches = 1000
inactive = 200
particles = 10000000


#########   Exporting to OpenMC geometry.xml and materials.xml Files  ##########

# Write all BEAVRS materials to materials.xml file
infermc.beavrs.make_iso_in_lab()
infermc.beavrs.write_materials_file()

# Extract reflected geometry from InferMC's pre-built assembly Geometries
full_core = infermc.beavrs.build_2d_core()
openmc_geometry = opencg_compatible.get_openmc_geometry(full_core)
openmc_geometry.export_to_xml()


##################   Exporting to OpenMC settings.xml File  ###################

# Construct uniform initial source distribution over fissionable zones
lower_left = [0., 0., full_core.bounds[2]]
upper_right = [+200., +200., full_core.bounds[5]]
source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
source.space.only_fissionable = True

settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.statepoint_interval = 50
settings_file.ptables = True
settings_file.output = {'tallies': False}
settings_file.source = source
settings_file.sourcepoint_write = False
settings_file.export_to_xml()


##################   Exporting to OpenMC plots.xml File  ######################

# Initialize the BEAVRS color mapping scheme
b = infermc.beavrs.beavrs
b.write_openmc_plots()

plot = openmc.Plot(plot_id=1)
bounds = full_core.bounds
plot.width = [full_core.max_x - full_core.min_x,
              full_core.max_y - full_core.min_y]
plot.origin = [bounds[0] + (bounds[3] - bounds[0]) / 2.,
               bounds[1] + (bounds[4] - bounds[1]) / 2.,
               bounds[2] + (bounds[5] - bounds[2]) / 2.]
plot.color = 'mat'
plot.filename = 'full-core'
plot.col_spec = b.plots.colspec_mat
plot.pixels = [5000, 5000]

plot_file = openmc.Plots([plot])
plot_file.export_to_xml()


########################   Build OpenMC MGXS Library  #########################

# Get all cells filled with a "fuel" material
mat_cells = openmc_geometry.get_all_material_cells()
fuel_cells = []
for cell in mat_cells:
    if 'fuel' in cell.fill.name.lower():
        fuel_cells.append(cell)

# Initialize a fine (40-) group "distribcell" MGXS Library for OpenMOC
cell_mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
cell_mgxs_lib.energy_groups = group_structures['CASMO']['40-group']
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
        mgxs.nuclides = ['U235', 'U238', 'total']

# Initialize a fine (40-) group "material" MGXS Library for OpenMOC
mat_mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
mat_mgxs_lib.energy_groups = group_structures['CASMO']['40-group']
mat_mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                           'chi', 'absorption', 'capture']
mat_mgxs_lib.domain_type = 'material'
mat_mgxs_lib.correction = None
mat_mgxs_lib.build_library()


###################  Create Mesh Tallies for Verification  ####################

# Instantiate a tally Mesh
mesh = openmc.Mesh(name='assembly mesh')
mesh.type = 'regular'
mesh.dimension = [int(np.ceil(15./2.*17)), int(np.ceil(15./2.*17)), 1]
mesh.lower_left = [-1.26492/2., -1.26492/2., 203.]
mesh.width = (1.26492, 1.26492, 10)

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
