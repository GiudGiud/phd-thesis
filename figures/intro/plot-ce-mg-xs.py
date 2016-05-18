import os
import re
import copy

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from infermc.energy_groups import group_structures
import openmc
import openmc.mgxs as mgxs
import pyne.ace

sns.set_style('ticks')

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
settings_file.particles = 100000

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

# Fine energy flux Tally
flux = openmc.Tally(name='flux')
energies = np.logspace(np.log10(1E-9), np.log10(10.), 2000)
flux.filters = [openmc.Filter(type='energy', bins=energies)]
flux.scores = ['flux']

# Add OpenMC tallies to the tallies file for XML generation
tallies_file = openmc.Tallies([flux])
tallies_file += list(fission.tallies.values())
tallies_file += list(capture.tallies.values())
tallies_file.export_to_xml()

# Run OpenMC
#openmc.run(output=True, mpi_procs=16, threads=1)

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')

# Load statepoint data into MGXS
fission.load_from_statepoint(sp)
capture.load_from_statepoint(sp)

# Get flux mean array and lengthen it to match energy bins for the plot
flux = sp.get_tally(name='flux')
flux_mean = flux.mean.flatten()
flux_mean = np.insert(flux_mean, 0, flux_mean[0])

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


for num_groups in [70, 16, 2]:

    # Condense to the coarse group structure
    energy_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
    fission = fission.get_condensed_xs(energy_groups)
    capture = capture.get_condensed_xs(energy_groups)

    ###########################  U-235 FISSION  ###############################
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plot1 = ax1.loglog(u235.energy, ce_fission.sigma, color='b',
                       linewidth=1, label='Continuous')

    # Extract energy group bounds and MGXS values to plot
    x = copy.deepcopy(fission.energy_groups.group_edges)
    y = fission.get_xs(
        nuclides=['U-235'], order_groups='decreasing',xs_type='micro')

    # Fix low energy bound to the value defined by the ACE library
    x[0] = u235.energy[0]
    y = np.insert(y, 0, y[0])

    # Create a step plot for the MGXS
    plot2 = ax1.plot(x, y, drawstyle='steps', color='r', linewidth=2,
                     label='{}-Group'.format(num_groups))

    # Begin plot customization
    ax1.set_xlabel('Energy [MeV]', fontsize=12)
    ax1.set_ylabel('Microscopic Capture XS [barns]', fontsize=12)

    # Plot the flux
    ax2 = ax1.twinx()
    plot3 = ax2.loglog(flux.filters[0].bins, flux_mean,
                       color='g', linewidth=1, zorder=1, label='Flux')
    ax2.set_ylabel('Flux', color='g', fontsize=12)
    ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
    ax1.patch.set_visible(False) # hide the 'canvas'

    # Create legend
    plots = plot1+plot2+plot3
    labels = [plot.get_label() for plot in plots]
    ax1.legend(plots, labels, loc='center right', fontsize=12)

    # Customize the plot
    sns.set_style('ticks')
    plt.title('U-235 Fission Cross Section', y=1.03, fontsize=16)
    plt.xlim((1e-9, 1e1))
    plt.savefig('u235-fission-{}.png'.format(num_groups), bbox_inches='tight')


    ###########################  U-238 CAPTURE  ###############################
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plot1 = ax1.loglog(u238.energy, ce_capture.sigma, color='b',
                       linewidth=1, label='Continuous')

    # Extract energy group bounds and MGXS values to plot
    x = copy.deepcopy(capture.energy_groups.group_edges)
    y = capture.get_xs(nuclides=['U-238'], order_groups='decreasing', xs_type='micro')

    # Fix low energy bound to the value defined by the ACE library
    x[0] = u238.energy[0]
    y = np.insert(y, 0, y[0])

    # Create a step plot for the MGXS
    plot2 = ax1.plot(x, y, drawstyle='steps', color='r',
                     linewidth=2, label='{}-Group'.format(num_groups))

    # Begin plot customization
    ax1.set_xlabel('Energy [MeV]', fontsize=12)
    ax1.set_ylabel('Microscopic Capture XS [barns]', fontsize=12)

    # Plot the flux
    ax2 = ax1.twinx()
    plot3 = ax2.loglog(flux.filters[0].bins, flux_mean,
               color='g', linewidth=1, zorder=1, label='Flux')
    ax2.set_ylabel('Flux', color='g', fontsize=12)
    ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
    ax1.patch.set_visible(False) # hide the 'canvas'

    # Create legend
    plots = plot1+plot2+plot3
    labels = [plot.get_label() for plot in plots]
    ax1.legend(plots, labels, loc='center right', fontsize=12)

    # Customize the plot
    sns.set_style('ticks')
    plt.title('U-238 Capture Cross Section', y=1.03, fontsize=16)
    plt.xlim((1e-9, 1e1))
    ax1.set_ylim((1E-4, 1e8))
    plt.savefig('u238-capture-{}.png'.format(num_groups), bbox_inches='tight')