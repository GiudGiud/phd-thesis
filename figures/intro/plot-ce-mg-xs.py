import os
import re

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from infermc.energy_groups import group_structures
import openmc
import openmc.mgxs as mgxs
import pyne.ace

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')
u238 = openmc.Nuclide('U-238')
zr90 = openmc.Nuclide('Zr-90')

# 1.6% enriched fuel
fuel = openmc.Material(name='1.6% Fuel')
fuel.set_density('g/cm3', 10.31341)
fuel.add_nuclide(u235, 3.7503e-4)
fuel.add_nuclide(u238, 2.2625e-2)
fuel.add_nuclide(o16, 4.6007e-2)

# borated water
water = openmc.Material(name='Borated Water')
water.set_density('g/cm3', 0.740582)
water.add_nuclide(h1, 4.9457e-2)
water.add_nuclide(o16, 2.4732e-2)

# zircaloy
zircaloy = openmc.Material(name='Zircaloy')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_nuclide(zr90, 7.2758e-3)

# Instantiate a Materials collection
materials_file = openmc.Materials((fuel, water, zircaloy))
materials_file.default_xs = '71c'

# Export to "materials.xml"
materials_file.export_to_xml()

# Create cylinders for the fuel and clad
fuel_outer_radius = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.39218)
clad_outer_radius = openmc.ZCylinder(x0=0.0, y0=0.0, R=0.45720)

# Create boundary planes to surround the geometry
min_x = openmc.XPlane(x0=-0.63, boundary_type='reflective')
max_x = openmc.XPlane(x0=+0.63, boundary_type='reflective')
min_y = openmc.YPlane(y0=-0.63, boundary_type='reflective')
max_y = openmc.YPlane(y0=+0.63, boundary_type='reflective')
min_z = openmc.ZPlane(z0=-0.63, boundary_type='reflective')
max_z = openmc.ZPlane(z0=+0.63, boundary_type='reflective')

# Create a Universe to encapsulate a fuel pin
pin_cell_universe = openmc.Universe(name='1.6% Fuel Pin')

# Create fuel Cell
fuel_cell = openmc.Cell(name='1.6% Fuel')
fuel_cell.fill = fuel
fuel_cell.region = -fuel_outer_radius
pin_cell_universe.add_cell(fuel_cell)

# Create a clad Cell
clad_cell = openmc.Cell(name='1.6% Clad')
clad_cell.fill = zircaloy
clad_cell.region = +fuel_outer_radius & -clad_outer_radius
pin_cell_universe.add_cell(clad_cell)

# Create a moderator Cell
moderator_cell = openmc.Cell(name='1.6% Moderator')
moderator_cell.fill = water
moderator_cell.region = +clad_outer_radius
pin_cell_universe.add_cell(moderator_cell)

# Create root Cell
root_cell = openmc.Cell(name='root cell')
root_cell.region = +min_x & -max_x & +min_y & -max_y
root_cell.fill = pin_cell_universe

# Create root Universe
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# Create Geometry and set root Universe
openmc_geometry = openmc.Geometry()
openmc_geometry.root_universe = root_universe
openmc_geometry.export_to_xml()

# Instantiate a Settings object
settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 10
settings_file.particles = 1000

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.63, -0.63, -0.63, 0.63, 0.63, 0.63]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

# Export to "settings.xml"
settings_file.export_to_xml()

# Extract all Cells filled by Materials
openmc_cells = openmc_geometry.get_all_material_cells()

# Create MGXS of interest
fission = mgxs.FissionXS(groups=group_structures['CASMO']['70-group'])
capture = mgxs.CaptureXS(groups=group_structures['CASMO']['70-group'])
fission.domain = fuel_cell
capture.domain = fuel_cell
fission.by_nuclide = True
capture.by_nuclide = True

# Add OpenMC tallies to the tallies file for XML generation
tallies_file = openmc.Tallies()
tallies_file += list(fission.tallies.values())
tallies_file += list(capture.tallies.values())
tallies_file.export_to_xml()

# Run OpenMC
openmc.run(output=True)

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')

# Load statepoint data into MGXS
fission.load_from_statepoint(sp)
capture.load_from_statepoint(sp)

# Instantiate a PyNE ACE continuous-energy cross sections library
xs_dir = os.environ['OPENMC_CROSS_SECTIONS']
xs_dir = re.sub('cross_sections.xml$', '', xs_dir)
u238_lib = pyne.ace.Library(xs_dir + '293.6K/U_238_293.6K.ace')
u235_lib = pyne.ace.Library(xs_dir + '293.6K/U_235_293.6K.ace')
u238_lib.read('92238.71c')
u235_lib.read('92235.71c')

# Extract the U-235 data from the library
u238 = u238_lib.tables['92238.71c']
u235 = u235_lib.tables['92235.71c']

# Extract the continuous-energy cross section data
ce_capture = u238.reactions[102]
ce_fission = u235.reactions[18]

# Create a loglog plot of the U-235 continuous-energy fission cross section
fig = plt.figure()
plt.loglog(u235.energy, ce_fission.sigma, color='b', linewidth=1)

# Extract energy group bounds and MGXS values to plot
x = fission.energy_groups.group_edges
y = fission.get_xs(nuclides=['U-235'], order_groups='decreasing', xs_type='micro')

# Fix low energy bound to the value defined by the ACE library
x[0] = u235.energy[0]
y = np.insert(y, 0, y[0])

# Create a step plot for the MGXS
plt.plot(x, y, drawstyle='steps', color='r', linewidth=3)
plt.title('U-235 Fission Cross Section')
plt.xlabel('Energy [MeV]')
plt.ylabel('Microscopic Fission XS [barns]')
plt.legend(['Continuous', '70-Group'])
plt.xlim((x.min(), x.max()))
plt.savefig('u235-fission-70.png', bbox_inches='tight')



# Create a loglog plot of the U-235 continuous-energy fission cross section
fig = plt.figure()
plt.loglog(u238.energy, ce_capture.sigma, color='b', linewidth=1)

# Extract energy group bounds and MGXS values to plot
x = capture.energy_groups.group_edges
y = capture.get_xs(nuclides=['U-238'], order_groups='decreasing', xs_type='micro')

# Fix low energy bound to the value defined by the ACE library
x[0] = u238.energy[0]
y = np.insert(y, 0, y[0])

# Create a step plot for the MGXS
plt.plot(x, y, drawstyle='steps', color='r', linewidth=3)
plt.title('U-238 Capture Cross Section')
plt.xlabel('Energy [MeV]')
plt.ylabel('Microscopic Capture XS [barns]')
plt.legend(['Continuous', '70-Group'])
plt.xlim((x.min(), x.max()))
plt.savefig('u238-capture-70.png', bbox_inches='tight')



# Condense to the 2-group structure
fission = fission.get_condensed_xs(group_structures['CASMO']['2-group'])
capture = capture.get_condensed_xs(group_structures['CASMO']['2-group'])


# Create a loglog plot of the U-235 continuous-energy fission cross section
fig = plt.figure()
plt.loglog(u235.energy, ce_fission.sigma, color='b', linewidth=1)

# Extract energy group bounds and MGXS values to plot
x = fission.energy_groups.group_edges
y = fission.get_xs(nuclides=['U-235'], order_groups='decreasing', xs_type='micro')

# Fix low energy bound to the value defined by the ACE library
x[0] = u235.energy[0]
y = np.insert(y, 0, y[0])

# Create a step plot for the MGXS
plt.plot(x, y, drawstyle='steps', color='r', linewidth=3)
plt.title('U-235 Fission Cross Section')
plt.xlabel('Energy [MeV]')
plt.ylabel('Microscopic Fission XS [barns]')
plt.legend(['Continuous', '2-Group'])
plt.xlim((x.min(), x.max()))
plt.savefig('u235-fission-2.png', bbox_inches='tight')



# Create a loglog plot of the U-235 continuous-energy capture cross section
fig = plt.figure()
plt.loglog(u238.energy, ce_capture.sigma, color='b', linewidth=1)

# Extract energy group bounds and MGXS values to plot
x = capture.energy_groups.group_edges
y = capture.get_xs(nuclides=['U-238'], order_groups='decreasing', xs_type='micro')

# Fix low energy bound to the value defined by the ACE library
x[0] = u238.energy[0]
y = np.insert(y, 0, y[0])

# Create a step plot for the MGXS
plt.plot(x, y, drawstyle='steps', color='r', linewidth=3)
plt.title('U-238 Capture Cross Section')
plt.xlabel('Energy [MeV]')
plt.ylabel('Microscopic Capture XS [barns]')
plt.legend(['Continuous', '16-Group'])
plt.xlim((x.min(), x.max()))
plt.savefig('u238-capture-16.png', bbox_inches='tight')
