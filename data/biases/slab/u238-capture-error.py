import numpy as np

import openmc.mgxs
import openmoc
from openmoc.opencg_compatible import get_openmoc_geometry
from infermc.energy_groups import group_structures


def get_fluxes(solver, mgxs_lib):
    """

    :param solver:
    :param mgxs_lib:
    :return:
    """

    # Extract the OpenMOC scalar fluxes
    openmoc_fluxes = openmoc.process.get_scalar_fluxes(solver)

    # Extract the OpenMC scalar fluxes
    openmoc_geometry = solver.getGeometry()
    num_fsrs = openmoc_geometry.getNumFSRs()
    num_groups = openmoc_geometry.getNumEnergyGroups()

    openmc_fluxes = np.zeros((num_fsrs, num_groups), dtype=np.float64)
    fsr_volumes = np.zeros(num_fsrs, dtype=np.float64)

    # Get the OpenMC flux in each FSR
    for fsr in range(num_fsrs):

        # Find the OpenMOC cell and volume for this FSR
        openmoc_cell = openmoc_geometry.findCellContainingFSR(fsr)
        cell_id = openmoc_cell.getId()
        fsr_volumes[fsr] = track_generator.getFSRVolume(fsr)

        # Store the volume-averaged flux
        mgxs = mgxs_lib.get_mgxs(cell_id, 'nu-fission')
        flux = mgxs.tallies['flux'].mean.flatten()
        flux = np.flipud(flux) / fsr_volumes[fsr]
        openmc_fluxes[fsr, :] = flux

    # Extract energy group edges
    group_edges = mgxs_lib.energy_groups.group_edges
    group_edges *= 1e6      # Convert to units of eV
    group_edges[0] = 1e-5     # Adjust lower bound (for loglog scaling)

    # Compute difference in energy bounds for each group
    group_deltas = np.ediff1d(group_edges)
    group_deltas = np.flipud(group_deltas)

    # Normalize fluxes to the total integrated flux
    openmc_fluxes /= np.sum(openmc_fluxes * group_deltas, axis=1)[:,np.newaxis]
    openmoc_fluxes /= np.sum(openmoc_fluxes * group_deltas, axis=1)[:,np.newaxis]
    return openmc_fluxes, openmoc_fluxes, fsr_volumes


openmoc.log.set_log_level('RESULT')
opts = openmoc.options.Options()

groups = [1, 2, 4] #, 8, 16, 25, 40, 70]
scatter = 'iso-in-lab'
num_mesh = 16

# Initialize a fine (70-)group MGXS Library from OpenMC statepoint data
directory = '{}/{}x/'.format(scatter, num_mesh)
sp = openmc.StatePoint(directory + 'statepoint.100.h5')
mgxs_lib = openmc.mgxs.Library.load_from_file(directory=directory)

abs_rel_err = np.zeros((len(groups), 3), dtype=np.float)
capt_rel_err = np.zeros((len(groups), 3), dtype=np.float)
u238_frac = np.zeros((len(groups), 3), dtype=np.float)

for i, num_groups in enumerate(groups):
    print('# groups = {}'.format(num_groups))

    # Build a coarse group Library from the fine (70-)group Library
    coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]
    condense_lib = mgxs_lib.get_condensed_library(coarse_groups)

    # Create an OpenMOC Geometry from the OpenCG Geometry
    openmoc_geometry = get_openmoc_geometry(condense_lib.opencg_geometry)
    openmoc.materialize.load_openmc_mgxs_lib(condense_lib, openmoc_geometry)

    # Generate tracks
    track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
    track_generator.setNumThreads(opts.num_omp_threads)
    track_generator.generateTracks()

    # Instantiate a Solver
    solver = openmoc.CPUSolver(track_generator)
    solver.setNumThreads(opts.num_omp_threads)
    solver.setConvergenceThreshold(1E-7)

    # Run OpenMOC
    solver.computeEigenvalue(opts.max_iters)

    # Extract the normalized fluxes and FSR volumes
    openmc_fluxes, openmoc_fluxes, fsr_volumes = get_fluxes(solver, condense_lib)

    # Extract the OpenMC scalar fluxes
    num_fsrs = openmoc_geometry.getNumFSRs()
    num_groups = openmoc_geometry.getNumEnergyGroups()

    # Allocate arrays for absorption/capture rates by material and energy group
    openmc_abs = np.zeros(num_groups, dtype=np.float)
    openmoc_abs = np.zeros(num_groups, dtype=np.float)
    openmc_capt = np.zeros(num_groups, dtype=np.float)
    openmoc_capt = np.zeros(num_groups, dtype=np.float)

    for fsr in range(num_fsrs):

        # Find the OpenMOC cell and volume for this FSR
        openmoc_cell = openmoc_geometry.findCellContainingFSR(fsr)
        cell_id = openmoc_cell.getId()

        # Get the capture cross section for this cell
        abs_mgxs = condense_lib.get_mgxs(cell_id, 'absorption')
        capt_mgxs = condense_lib.get_mgxs(cell_id, 'capture')

        # Compute OpenMC/OpenMOC total capture rates
        abs_mean = abs_mgxs.get_xs(nuclides='sum', xs_type='macro')
        openmc_abs += openmc_fluxes[fsr,:] * abs_mean.flatten() * fsr_volumes[fsr]
        openmoc_abs += openmoc_fluxes[fsr,:] * abs_mean.flatten() * fsr_volumes[fsr]

        # Compute OpenMC/OpenMOC U-238 capture rates
        if openmoc_cell.getName() == 'fuel':
            capt_mean = capt_mgxs.get_xs(nuclides=['U-238'], xs_type='macro')
            openmc_capt += openmc_fluxes[fsr,:] * capt_mean.flatten() * fsr_volumes[fsr]
            openmoc_capt += openmoc_fluxes[fsr,:] * capt_mean.flatten() * fsr_volumes[fsr]

    # Compute the percent rel. err. in group 27
    abs_rel_err[i,0] = (openmoc_abs[26] - openmc_abs[26]) / openmc_abs[26] * 100.
    capt_rel_err[i,0] = (openmoc_capt[26] - openmc_capt[26]) / openmc_capt[26] * 100

    # Compute the percent rel. err. in groups 14-27
    abs_rel_err[i,1] = (np.sum(openmoc_abs[13:26]) - np.sum(openmc_abs[13:26])) / np.sum(openmc_abs[13:26]) * 100.
    capt_rel_err[i,1] = (np.sum(openmoc_capt[13:26]) - np.sum(openmc_capt[13:26])) / np.sum(openmc_capt[13:26]) * 100

    # Compute the percent rel. err. in all groups
    abs_rel_err[i,2] = (np.sum(openmoc_abs) - np.sum(openmc_abs)) / np.sum(openmc_abs) * 100.
    capt_rel_err[i,2] = (np.sum(openmoc_capt) - np.sum(openmc_capt)) / np.sum(openmc_capt) * 100

    # Compute the percentage of U-238 capture to total absorption
    u238_frac[i,0] = openmoc_capt[26] / openmoc_abs[26] * 100.
    u238_frac[i,1] = np.sum(openmoc_capt[13:26]) / np.sum(openmoc_abs[13:26]) * 100.
    u238_frac[i,2] = np.sum(openmoc_capt) / np.sum(openmoc_abs) * 100.

group_ranges = ['Group 27', 'Groups 14-27', 'All Groups']

# Print U-238 capture and total absorption table for LaTeX
print('Capture and Absorption Rel. Err.')
for i, num_groups in enumerate(groups):
    row = '{} &'.format(num_groups)
    for j, group_range in enumerate(group_ranges):
        row += ' {:1.2f} &'.format(capt_rel_err[i,j])
    row += ' &'
    for j, group_range in enumerate(group_ranges):
        row += ' {:1.2f} &'.format(abs_rel_err[i,j])
    print(row[:-1] + '\\\\')

# Print U-238 capture-to-total absorption table for LaTeX
print('Capture and Absorption Rel. Err.')
for i, num_groups in enumerate(groups):
    row = '{} &'.format(num_groups)
    for j, group_range in enumerate(group_ranges):
        row += ' {:2.2f} &'.format(u238_frac[i,j])
    print(row[:-1] + '\\\\')
