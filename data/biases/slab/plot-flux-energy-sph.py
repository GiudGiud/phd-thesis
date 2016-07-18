import os
import re

import numpy as np
import matplotlib

# Force headless backend for plotting on clusters
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns

import openmc.mgxs
import openmc.opencg_compatible
import openmoc
import openmoc.opencg_compatible
from infermc.energy_groups import group_structures
import pyne.ace

# Force non-interactive mode for plotting on clusters
plt.ioff()
sns.set_style('ticks')

openmoc.log.set_log_level('NORMAL')
opts = openmoc.options.Options()

# Query the user for the number of energy groups
scatter = 'iso-in-lab'
mesh = 1
num_groups = 70

directory = '{}/{}x'.format(scatter, mesh)

# Create an OpenMOC Geometry from the OpenCG Geometry
mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)
openmoc_geometry = \
    openmoc.opencg_compatible.get_openmoc_geometry(mgxs_lib.opencg_geometry)

coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
mgxs_lib = mgxs_lib.get_condensed_library(coarse_groups)

# FIXME
# Compute SPH factors
sph, sph_mgxs_lib, sph_indices = \
    openmoc.materialize.compute_sph_factors(
        mgxs_lib, num_azim=128, azim_spacing=0.01, sph_tol=1E-7,
        num_threads=opts.num_omp_threads, max_sph_iters=20)

# Load the SPH-corrected MGXS library data
materials = \
    openmoc.materialize.load_openmc_mgxs_lib(sph_mgxs_lib, openmoc_geometry)

# Initialize an OpenMOC TrackGenerator and Solver
track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.01)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()

# Initialize an OpenMOC Solver
solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(1E-7)
solver.setNumThreads(opts.num_omp_threads)

# Run an eigenvalue calulation with the MGXS from OpenMC
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()


###############################################################################
#                         Extracting Scalar Fluxes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

# Extract the OpenMOC scalar fluxes
openmoc_fluxes = openmoc.process.get_scalar_fluxes(solver) # * sph

# Extract the OpenMC scalar fluxes
num_fsrs = openmoc_geometry.getNumFSRs()
num_groups = openmoc_geometry.getNumEnergyGroups()
openmc_fluxes = np.zeros((num_fsrs, num_groups), dtype=np.float64)
fsr_volumes = np.zeros(num_fsrs, dtype=np.float64)
openmc_fiss = 0.

# Get innermost/outermost fuel FSRs
fuel_centroids = []
fuel_fsrs = []

# Get the OpenMC flux in each FSR
for fsr in range(num_fsrs):

    # Find the OpenMOC cell and volume for this FSR
    openmoc_cell = openmoc_geometry.findCellContainingFSR(fsr)
    cell_id = openmoc_cell.getId()
    fsr_volumes[fsr] = track_generator.getFSRVolume(fsr)

    # Update min/max fuel FSR
    if openmoc_cell.getName() == 'fuel':
        centroid = openmoc_geometry.getFSRCentroid(fsr)
        x = centroid.getX()
        fuel_centroids.append(x)
        fuel_fsrs.append(fsr)

    # Store the volume-averaged flux
    mgxs = mgxs_lib.get_mgxs(cell_id, 'nu-fission')
    mgxs_mean = mgxs.get_xs(nuclides='sum', xs_type='macro')
    flux = mgxs.tallies['flux'].mean.flatten()
    openmc_fluxes[fsr, :] = flux[::-1] / fsr_volumes[fsr]
    openmc_fiss += np.sum(flux[::-1] * mgxs_mean)

openmc_fluxes /= openmc_fiss

fuel_centroids = np.asarray(fuel_centroids)
fuel_fsrs = np.asarray(fuel_fsrs)
min_x = np.argmin(fuel_centroids)
max_x = np.argmin(fuel_centroids)

# Extract energy group edges
group_edges = mgxs_lib.energy_groups.group_edges
group_edges *= 1e6      # Convert to units of eV
group_edges[0] = 1e-5   # Adjust lower bound (for loglog scaling)

# FIXME
openmc_fluxes = openmc_fluxes[:, ::-1]
openmoc_fluxes = openmoc_fluxes[:, ::-1]

# Extend the mgxs values array for matplotlib's step plot of fluxes
openmc_fluxes = np.insert(openmc_fluxes, 0, openmc_fluxes[:,0], axis=1)
openmoc_fluxes = np.insert(openmoc_fluxes, 0, openmoc_fluxes[:,0], axis=1)


###############################################################################
#                 Plot OpenMC-to-OpenMOC Scalar Flux Errors
###############################################################################

# Compute volume-weighted FSR fluxes across the geometry
vol_avg_openmc_fluxes = np.nansum(openmc_fluxes[fuel_fsrs,:] *
                                  fsr_volumes[fuel_fsrs, np.newaxis], axis=0)
vol_avg_openmoc_fluxes = np.nansum(openmoc_fluxes[fuel_fsrs,:] *
                                   fsr_volumes[fuel_fsrs, np.newaxis], axis=0)

# Normalize the fluxes to the mean
vol_avg_openmc_fluxes /= np.mean(vol_avg_openmc_fluxes)
vol_avg_openmoc_fluxes /= np.mean(vol_avg_openmoc_fluxes)

# Compute lethargy for each group
lethargies = np.ediff1d(np.log(group_edges))

# Normalize the fluxes to the mean
vol_avg_openmc_fluxes[1:] /= lethargies
vol_avg_openmoc_fluxes[1:] /= lethargies

vol_avg_openmc_fluxes[0] = vol_avg_openmc_fluxes[1]
vol_avg_openmoc_fluxes[0] = vol_avg_openmoc_fluxes[1]

plt.plot(group_edges, vol_avg_openmc_fluxes,
         drawstyle='steps', color='r', linewidth=2)
plt.plot(group_edges, vol_avg_openmoc_fluxes,
         drawstyle='steps', color='b', linewidth=2)

plt.yscale('log')
plt.xscale('log')
plt.xlabel('Energy [eV]', fontsize=12)
plt.ylabel('Flux / Lethargy', fontsize=12)
plt.xlim((min(group_edges), max(group_edges)))
plt.legend(['OpenMC', 'OpenMOC'], loc='best', fontsize=12)
plt.savefig('vol-avg-flux-sph.png', bbox_inches='tight')
plt.close()

###############################################################################
#                 Plot OpenMC-to-OpenMOC Scalar Flux Errors
###############################################################################

# Instantiate a PyNE ACE continuous-energy cross sections library
xs_dir = os.environ['OPENMC_CROSS_SECTIONS']
xs_dir = re.sub('cross_sections.xml$', '', xs_dir)
pyne_lib = pyne.ace.Library(xs_dir + '293.6K/U_238_293.6K.ace')
pyne_lib.read('92238.71c')

# Extract the U-235 data from the library
u238 = pyne_lib.tables['92238.71c']

# Extract the continuous-energy U-238 fission cross section data
capture = u238.reactions[102]

# Compute the percent relative error in the flux
delta_flux = openmoc_fluxes - openmc_fluxes
rel_err = delta_flux / openmc_fluxes * 100.

'''
# Compute the percent relative error in the flux
rel_err = np.zeros(openmc_fluxes.shape)
for fsr in range(num_fsrs):
    delta_flux = openmoc_fluxes[fsr,:] - openmc_fluxes[fsr,:]
    rel_err[fsr,:] = delta_flux / openmc_fluxes[fsr,:] * 100.
'''

'''
# Plot OpenMOC relative flux errors in each FSR
for fsr in range(num_fsrs):

    # Get the OpenMOC cell and material for this FSR
    cell = openmoc_geometry.findCellContainingFSR(fsr)

    # Create a step plot for the MGXS
    fig, ax1 = plt.subplots()

    plt.plot(group_edges, rel_err[fsr,:],
             drawstyle='steps', color='r', linewidth=2)

    plt.xlabel('Energy [eV]', fontsize=12)
    plt.ylabel('Relative Error [%]', fontsize=12)
    plt.xlim((min(group_edges), max(group_edges)))
    plt.xscale('log')

    # Create loglog plot of U-238 continuous-energy capture cross section 
    if 'fuel' in cell.getName():
        ax2 = ax1.twinx()
        ax2.loglog(u238.energy*1e6, capture.sigma, color='g', \
                   linewidth=1, zorder=1)
        ax2.set_ylabel('U-238 Capture XS [barns]', color='g', fontsize=12)
        ax2.set_yscale('log')

        ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
        ax1.patch.set_visible(False) # hide the 'canvas' 

    filename = 'rel-err-{0}.png'.format(fsr)
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
'''

###############################################################################
#        Plot OpenMC-to-OpenMOC Innermost/Outermost Scalar Flux Error
###############################################################################

min_fsr = fuel_fsrs[np.argmin(fuel_centroids)]
max_fsr = fuel_fsrs[np.argmax(fuel_centroids)]

# Plot the error for the innermost and outermost FSRS atop each other
fig, ax1 = plt.subplots()

plt.plot(group_edges, rel_err[min_fsr,:],
         drawstyle='steps', color='b', linewidth=2)
plt.plot(group_edges, rel_err[max_fsr,:],
         drawstyle='steps', color='r', linewidth=2)

plt.xlabel('Energy [eV]', fontsize=12)
plt.ylabel('Relative Error [%]', fontsize=12)
plt.xlim((min(group_edges), max(group_edges)))
plt.xscale('log')
plt.legend(['Innermost', 'Outermost'], fontsize=12)

# Create loglog plot of U-238 continuous-energy capture cross section
ax2 = ax1.twinx()
ax2.loglog(u238.energy*1e6, capture.sigma, color='g', \
           linewidth=1, zorder=1)
ax2.set_ylabel('U-238 Capture XS [barns]', color='g', fontsize=12)
ax2.set_yscale('log')

ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2
ax1.patch.set_visible(False) # hide the 'canvas'

plt.savefig('rel-err-inner-outer-sph.png', bbox_inches='tight')
plt.close()


###############################################################################
#        Plot OpenMC-to-OpenMOC Flux Error Across Fuel FSRs in Group 27
###############################################################################

# Plot the relative error for group 27
fig = plt.figure()

# FIXME
group_index = 70 - 27 + 1

centroid_indices = np.argsort(fuel_centroids)
fuel_indices = fuel_fsrs[centroid_indices]
rel_err = rel_err[fuel_indices, group_index]
rel_err = np.insert(rel_err, 0, rel_err[0], axis=0)

plt.plot(np.arange(rel_err.shape[0]), rel_err, drawstyle='steps', color='b', linewidth=3)
plt.xlabel('Fuel FSR', fontsize=12)
plt.ylabel('Relative Error [%]', fontsize=12)
plt.xlim((0, rel_err.shape[0]-1))
plt.savefig('rel-err-fuel-fsrs-sph.png', bbox_inches='tight')
plt.close()


'''
###############################################################################
#        Plot OpenMC-to-OpenMOC Flux Error Across Fuel FSRs in Ranges A, B, C
###############################################################################

# Find energy group which encompasses 6.67 eV resonance
min_ind = mgxs_lib.energy_groups.get_group(6.67e-6) - 1

# Find energy group which encompasses 200 kEV
max_ind = mgxs_lib.energy_groups.get_group(2.e-2) - 1

# Compute OpenMOC flux in ranges A, B, and C for each fuel FSR
openmoc_flux_a = openmoc_fluxes[fuel_indices, group_index]
openmoc_flux_b = np.sum(openmoc_fluxes[fuel_indices, max_ind:min_ind], axis=1)
openmoc_flux_c = np.sum(openmoc_fluxes[fuel_indices, :], axis=1)

# Compute OpenMC flux in ranges A, B, and C for each fuel FSR
openmc_flux_a = openmc_fluxes[fuel_indices, group_index]
openmc_flux_b = np.sum(openmc_fluxes[fuel_indices, max_ind:min_ind], axis=1)
openmc_flux_c = np.sum(openmc_fluxes[fuel_indices, :], axis=1)

# Compute flux rel. err. for ranges A, B, and C for each fuel FSR
delta_flux_a = openmoc_flux_a - openmc_flux_a
delta_flux_b = openmoc_flux_b - openmc_flux_b
delta_flux_c = openmoc_flux_c - openmc_flux_c
rel_err_a = delta_flux_a / openmc_flux_a * 100.
rel_err_b = delta_flux_b / openmc_flux_b * 100.
rel_err_c = delta_flux_c / openmc_flux_c * 100.

# Plot the relative errors
fig = plt.figure()
plt.plot(np.arange(rel_err_a.shape[0]), rel_err_a, drawstyle='steps', color='r', linewidth=3)
plt.plot(np.arange(rel_err_b.shape[0]), rel_err_b, drawstyle='steps', color='g', linewidth=3)
plt.plot(np.arange(rel_err_c.shape[0]), rel_err_c, drawstyle='steps', color='b', linewidth=3)
plt.legend(['Range A', 'Range B', 'Range C'])
plt.xlabel('Fuel FSR', fontsize=12)
plt.ylabel('Relative Error [%]', fontsize=12)
plt.xlim((0, rel_err_a.shape[0]-1))
plt.savefig('rel-err-fuel-fsrs-ranges-sph.png', bbox_inches='tight')
plt.close()
'''