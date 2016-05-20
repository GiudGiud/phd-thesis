import os
import glob

import opencg
import openmc
import openmc.mgxs
import openmc.opencg_compatible
from infermc.energy_groups import group_structures

##################   Exporting to OpenMC materials.xml File  ###################

# Instantiate a Material and register the Nuclides
medium = openmc.Material(name='moderator')
medium.set_density('g/cc', 5.)
medium.add_nuclide('H-1',  0.028999667)
medium.add_nuclide('O-16', 0.01450188)
medium.add_nuclide('U-235', 0.000114142)
medium.add_nuclide('U-238', 0.006886019)
medium.add_nuclide('Zr-90', 0.002116053)

# Instantiate a MaterialsFile, register Material, and export to XML
materials_file = openmc.Materials([medium])
materials_file.default_xs = '71c'

# Get OpenCG versions of the OpenMC Material to fill OpenCG Cells below
opencg_medium = openmc.opencg_compatible.get_opencg_material(medium)


###################   Exporting to OpenMC geometry.xml File  ###################

# Instantiate boundary Planes
min_x = opencg.XPlane(boundary='reflective', x0=-0.63)
max_x = opencg.XPlane(boundary='reflective', x0=0.63)
min_y = opencg.YPlane(boundary='reflective', y0=-0.63)
max_y = opencg.YPlane(boundary='reflective', y0=0.63)

# Instantiate a Cell
root_cell = opencg.Cell(cell_id=1, name='cell')

# Register bounding Surfaces with the Cell
root_cell.add_surface(surface=min_x, halfspace=+1)
root_cell.add_surface(surface=max_x, halfspace=-1)
root_cell.add_surface(surface=min_y, halfspace=+1)
root_cell.add_surface(surface=max_y, halfspace=-1)

# Fill the Cell with the Material
root_cell.fill = opencg_medium

# Instantiate Universe
root_universe = opencg.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# Create OpenCG Geometry and set root Universe
opencg_geometry = opencg.Geometry()
opencg_geometry.root_universe = root_universe
openmc_geometry = openmc.opencg_compatible.get_openmc_geometry(opencg_geometry)


###################   Exporting to OpenMC settings.xml File  ###################

# Construct uniform initial source distribution over fissionable zones
lower_left = opencg_geometry.bounds[:3]
upper_right = opencg_geometry.bounds[3:]
source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
source.space.only_fissionable = True

# Create settings.xml files
settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 1
settings_file.particles = 100000000
settings_file.output = {'tallies': False}
settings_file.source = source
settings_file.sourcepoint_write = False


########################   Build OpenMC MGXS Library  #########################

# Initialize a fine (70-) group MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=True)
mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'cell'
mgxs_lib.correction = None
mgxs_lib.build_library()

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies_file, merge=True)


######################   Move Files into Directories  #########################

# Create files for anisotropic scattering
materials_file.export_to_xml()
openmc_geometry.export_to_xml()
settings_file.export_to_xml()
tallies_file.export_to_xml()

# Move files
for xml_file in glob.glob('*.xml'):
    os.rename(xml_file, 'anisotropic/{}'.format(xml_file))

# Create files for isotropic-in-lab
materials_file.make_isotropic_in_lab()
materials_file.export_to_xml()
openmc_geometry.export_to_xml()
settings_file.export_to_xml()
tallies_file.export_to_xml()

# Move files
for xml_file in glob.glob('*.xml'):
    os.rename(xml_file, 'iso-in-lab/{}'.format(xml_file))
