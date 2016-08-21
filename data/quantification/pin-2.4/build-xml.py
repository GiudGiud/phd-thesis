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

# Extract fuel pin of interest from InferMC's pre-built pin cell Geometries
pin_geometry = infermc.beavrs.find_pin('Fuel rod active region - 2.4% enr')
openmc_geometry = opencg_compatible.get_openmc_geometry(pin_geometry)
openmc_geometry.export_to_xml()


##################   Exporting to OpenMC settings.xml File  ###################

# Construct uniform initial source distribution over fissionable zones
lower_left = pin_geometry.bounds[:3]
upper_right = pin_geometry.bounds[3:]
lower_left[-1] = -10.
upper_right[-1] = 10.
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


#########################   Build OpenMC MGXS Library  #########################

# Initialize a fine (70-) group MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mgxs_lib.mgxs_types = ['total', 'fission', 'nu-fission', 'nu-scatter matrix',
                       'chi', 'absorption', 'capture']
mgxs_lib.domain_type = 'material'
mgxs_lib.correction = None
mgxs_lib.build_library()

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()
