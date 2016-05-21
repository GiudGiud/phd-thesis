import os
import glob
import copy

import opencg
import openmc
import openmc.mgxs
import openmc.opencg_compatible
from infermc.energy_groups import group_structures

##################   Exporting to OpenMC materials.xml File  ###################

# 3.1% w/o enriched UO2
fuel = openmc.Material(name='3.1 w/o enriched UO2')
fuel.set_density('sum')
fuel.add_nuclide('U-235', 7.18132E-4,)
fuel.add_nuclide('U-238', 2.21546E-2)
fuel.add_nuclide('O-16', 4.57642E-2)

# Cladding
clad = openmc.Material(name='Cladding')
clad.set_density('sum')
clad.add_nuclide('Zr-90', 2.18865E-2)
clad.add_nuclide('Zr-91', 4.77292E-3)
clad.add_nuclide('Zr-92', 7.29551E-3)
clad.add_nuclide('Zr-94', 7.39335E-3)
clad.add_nuclide('Zr-96', 1.19110E-3)

# Water
water = openmc.Material(name='Borated Water')
water.set_density('sum')
water.add_nuclide('H-1', 4.41459E-2)
water.add_nuclide('O-16', 2.20729E-2)
water.add_nuclide('B-10', 9.52537E-6)
water.add_nuclide('B-11', 3.83408E-5)
water.add_s_alpha_beta(name='HH2O', xs='71t')

# Instantiate a MaterialsFile, add Materials
materials_file = openmc.Materials((fuel, clad, water))
materials_file.default_xs = '71c'

# Get OpenCG versions of each OpenMC material to fill OpenCG Cells below
opencg_fuel = openmc.opencg_compatible.get_opencg_material(fuel)
opencg_clad = openmc.opencg_compatible.get_opencg_material(clad)
opencg_water = openmc.opencg_compatible.get_opencg_material(water)


###################   Exporting to OpenMC geometry.xml File  ###################

# Create bounding surfaces
min_x = opencg.XPlane(boundary='reflective', x0=0.0)
max_x = opencg.XPlane(boundary='reflective', x0=5.0)
min_y = opencg.YPlane(boundary='reflective', y0=0.0)
max_y = opencg.YPlane(boundary='reflective', y0=10.0)
min_z = opencg.ZPlane(boundary='reflective', z0=0.0)
max_z = opencg.ZPlane(boundary='reflective', z0=10.0)

# Create material interfacial surfaces
left = opencg.XPlane(surface_id=1, boundary='interface', x0=2.0)
right = opencg.XPlane(surface_id=2, boundary='interface', x0=2.4)

# Create a Universe to encapsulate the 1D slab
slab_universe = opencg.Universe(name='1D slab')

# Create fuel Cell
fuel_cell = opencg.Cell(name='fuel')
fuel_cell.fill = opencg_fuel
fuel_cell.add_surface(halfspace=+1, surface=min_x)
fuel_cell.add_surface(halfspace=-1, surface=left)
slab_universe.add_cell(fuel_cell)

# Create clad Cell
clad_cell = opencg.Cell(name='clad')
clad_cell.fill = opencg_clad
clad_cell.add_surface(halfspace=+1, surface=left)
clad_cell.add_surface(halfspace=-1, surface=right)
slab_universe.add_cell(clad_cell)

# Create water Cell
water_cell = opencg.Cell(name='water')
water_cell.fill = opencg_water
water_cell.add_surface(halfspace=+1, surface=right)
water_cell.add_surface(halfspace=-1, surface=max_x)
slab_universe.add_cell(water_cell)

# Create root Cell
root_cell = opencg.Cell(name='root cell')
root_cell.fill = slab_universe

# Add boundary planes
root_cell.add_surface(halfspace=+1, surface=min_x)
root_cell.add_surface(halfspace=-1, surface=max_x)
root_cell.add_surface(halfspace=+1, surface=min_y)
root_cell.add_surface(halfspace=-1, surface=max_y)
root_cell.add_surface(halfspace=+1, surface=min_z)
root_cell.add_surface(halfspace=-1, surface=max_z)

# Create root Universe
root_universe = opencg.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# Create Geometry and set root Universe
opencg_geometry = opencg.Geometry()
opencg_geometry.root_universe = root_universe

# Get an OpenMC version of this OpenCG geometry
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
settings_file.inactive = 10
settings_file.particles = 10000000
settings_file.output = {'tallies': False}
settings_file.source = source
settings_file.sourcepoint_write = False


####################   Exporting to OpenMC plots.xml File  #####################

# Find the edge widths for the geometry's bounding box
delta_x = opencg_geometry.max_x - opencg_geometry.min_x
delta_y = opencg_geometry.max_y - opencg_geometry.min_y
delta_z = opencg_geometry.max_z - opencg_geometry.min_z

# Instantiate a Plot
plot = openmc.Plot()
plot.origin = [delta_x/2., delta_y/2., delta_z/2.]
plot.width = [delta_x, delta_y]
plot.pixels = [1000, 250]
#plot.col_spec = {fuel.id: (255, 0, 0),
#                 clad.id: (120, 120, 120),
#                 water.id: (0, 0, 255)}
#plot.color = 'mat'
plot.color = 'cell'

# Instantiate a PlotsFile, add Plot, and export to "plots.xml"
plots_file = openmc.Plots([plot])


######################   Move Files into Directories  #########################

scattering = ['anisotropic', 'transport', 'iso-in-lab']
mesh = [1, 2, 4, 8, 16, 32, 64]

for scatter in scattering:
    print(scatter)
    for num_mesh in mesh:
        print('# mesh: {}'.format(num_mesh))

        if scatter == 'iso-in-lab':
            materials_file.make_isotropic_in_lab()
        materials_file.export_to_xml()

        # Copy the slab universe for meshing
        copy_slab_univ = copy.deepcopy(slab_universe)
        root_cell.fill = copy_slab_univ

        all_cells = copy_slab_univ.get_all_cells()
        for cell_id, cell in all_cells.items():
            linear_mesh = opencg.LinearMesh('x', cell.min_x, cell.max_x, num_mesh)
            new_cells = linear_mesh.subdivide_cell(cell, copy_slab_univ)

        # Get an OpenMC version of this OpenCG geometry
        openmc_geometry = openmc.opencg_compatible.get_openmc_geometry(opencg_geometry)
        openmc_geometry.export_to_xml()

        settings_file.export_to_xml()
        plots_file.export_to_xml()

        # Initialize a fine (70-) group MGXS Library for OpenMOC
        mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=False)
        mgxs_lib.energy_groups = group_structures['CASMO']['70-group']
        if scatter == 'transport':
            mgxs_lib.mgxs_types = ['nu-transport', 'nu-fission', 'nu-scatter matrix', 'chi']
            mgxs_lib.correction = 'P0'
        else:
            mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
            mgxs_lib.correction = None
        mgxs_lib.domain_type = 'cell'
        mgxs_lib.build_library()

        # Create a "tallies.xml" file for the MGXS Library
        tallies_file = openmc.Tallies()
        mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
        tallies_file.export_to_xml()

        # Move files
        for xml_file in glob.glob('*.xml'):
            os.rename(xml_file, '{}/{}x/{}'.format(scatter, num_mesh, xml_file))