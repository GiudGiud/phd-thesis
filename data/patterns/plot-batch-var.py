import os
import sys
import glob

import matplotlib

# force headless backend, or set 'backend' to 'Agg'
#  in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from _collections import OrderedDict

import openmc
import openmc.mgxs
import infermc
from infermc.energy_groups import group_structures


sns.set_style('ticks')


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directory = sys.argv[1]
nuclide = sys.argv[2]
num_groups = int(sys.argv[3])
sp_start = int(sys.argv[4])
sp_stride = int(sys.argv[5])

print('Directory: {}'.format(directory))
print('Nuclide: {}'.format(nuclide))
print('# Groups: {}'.format(num_groups))

if directory == 'full-core':
    directories = OrderedDict({directory: 'Full Core'})
else:
    directories = OrderedDict({directory: '1.6% Enr. (no BPs)'})

coarse_groups = group_structures['CASMO']['{}-group'.format(num_groups)]

if nuclide == 'U238':
    mgxs_type = 'capture'
else:
    mgxs_type = 'fission'

cwd = os.getcwd()

for directory in directories:
    print(directory)

    os.chdir(os.path.join(cwd, directory))

    sp = openmc.StatePoint('statepoint.10000.h5')
    statepoints = sorted(glob.glob('statepoint.*.h5'))[sp_start::sp_stride]

    if 'statepoint.10000.h5' not in statepoints:
        statepoints.append('statepoint.10000.h5')

    # Load the MGXS library for this benchmark
    mgxslib = openmc.mgxs.Library.load_from_file(
        directory='mgxs', filename='distribcell')
    mgxslib.load_from_statepoint(sp)

    if nuclide == 'U238':
        group = coarse_groups.get_group(6.67e-6)
    else:
        group = coarse_groups.get_group(1.e-11)
    lower_bound = coarse_groups.get_group_bounds(group)[0]

    # Loop over all fuel domains for each plot type
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))

        mgxs = mgxslib.get_mgxs(domain, mgxs_type)

        # Create cluster IDs from OpenMC subdomains
        clusters = infermc.featurize.get_lns(
            mgxs, mgxslib.openmc_geometry, repeat=False)

        # FIXME: Collect LNS together
        subdomains = []
        for cluster in np.unique(clusters):
            where = np.where(clusters == cluster)[0]
            subdomains.append(where)

        lns_mgxs_pf = infermc.series.get_mgxs_batch_series(
            mgxs, subdomains=subdomains, coarse_groups=coarse_groups,
            xs_type='micro', nuclides=[nuclide], groups=[group],
            statepoints=statepoints)

        print('LNS DONE')

        degen_mgxs_pf = infermc.series.get_mgxs_batch_series(
            mgxs, coarse_groups=coarse_groups, xs_type='micro',
            nuclides=[nuclide], groups=[group],
            statepoints=statepoints)

        print('DEGEN DONE')

        null_mgxs_pf = infermc.series.get_mgxs_batch_series(
            mgxs, subdomains=['all'], coarse_groups=coarse_groups,
            xs_type='micro', nuclides=[nuclide], groups=[group],
            statepoints=statepoints)

        print('NULL DONE')

        degen_df = degen_mgxs_pf.minor_xs('std. dev.')
        degen_df = degen_df.transpose()
        degen_df.index -= sp.n_inactive
        degen_df.index *= 1e-1

        null_df = null_mgxs_pf.minor_xs('std. dev.')
        null_df = null_df.transpose()
        null_df.index -= sp.n_inactive
        null_df.index *= 1e-1

        lns_df = lns_mgxs_pf.minor_xs('std. dev.')
        lns_df = lns_df.transpose()
        lns_df.index -= sp.n_inactive
        lns_df.index *= 1e-1
        
        # DEGENERATE
        degen_df = degen_df.abs()
        mean_degen_df = degen_df.mean(axis=1)
        max_degen_df = degen_df.max(axis=1)

        # NULL
        null_df = null_df.abs()

        # LNS
        lns_df = lns_df.abs()
        mean_lns_df = lns_df.mean(axis=1)
        max_lns_df = lns_df.max(axis=1)

        # Make filename
        if directory in ['2x2', 'reflector', 'full-core']:
            if '1.6' in domain.fill.name:
                filename = '16-enr-var-{}-{}.png'.format(mgxs_type, group)
            elif '2.4' in domain.fill.name:
                filename = '24-enr-var-{}-{}.png'.format(mgxs_type, group)
            elif '3.1' in domain.fill.name:
                filename = '31-enr-var-{}-{}.png'.format(mgxs_type, group)
        else:
            filename = '{}-var-{}-{}.png'.format(directory.replace('.', ''), mgxs_type, group)

        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'convergence/')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass

        # Plot the figure
        plt.figure()
        null_df.plot(loglog=True, legend=False, linestyle='-', color='g')
        max_lns_df.plot(loglog=True, legend=False, linestyle='-', color='r')
        mean_lns_df.plot(loglog=True, legend=False, linestyle=':', color='r')
        max_degen_df.plot(loglog=True, legend=False, linestyle='-', color='b')
        mean_degen_df.plot(loglog=True, legend=False, linestyle=':', color='b')

        # FIXME: Line for NULL
        N = sp.n_batches - sp.n_inactive - 20
        x0 = null_df.index.values[0]
        x1 = null_df.index.values[-1]
        y0 = null_df.iloc[0].values[0]
        N = np.float(x1) / np.float(x0)
        y1 = 1. / np.sqrt(N) * y0
        print('NULL: (x0,y0) = ({},{}), (x1,y1) = ({},{})'.format(x0,y0,x1,y1))
        plt.plot([x0, x1], [y0, y1], linestyle='-.', color='k', linewidth=1)

        # FIXME: Line for DEGENERATE (MAX)
        y0 = max_degen_df.iloc[0]
        y1 = 1. / np.sqrt(N) * y0
        print('DEGENERATE (MAX): (x0,y0) = ({},{}), (x1,y1) = ({},{})'.format(x0,y0,x1,y1))
        plt.plot([x0, x1], [y0, y1], linestyle='-.', color='k', linewidth=1)

        # FIXME: Line for DEGENERATE (MEAN)
        y0 = mean_degen_df.iloc[0]
        y1 = 1. / np.sqrt(N) * y0
        print('DEGENERATE (MEAN): (x0,y0) = ({},{}), (x1,y1) = ({},{})'.format(x0,y0,x1,y1))

        # FIXME: Line for LNS (MAX)
        y0 = max_lns_df.iloc[0]
        y1 = 1. / np.sqrt(N) * y0
        print('LNS (MAX): (x0,y0) = ({},{}), (x1,y1) = ({},{})'.format(x0,y0,x1,y1))
        plt.plot([x0, x1], [y0, y1], linestyle='-.', color='k', linewidth=1)

        # FIXME: Line for LNS (MEAN)
        y0 = mean_lns_df.iloc[0]
        y1 = 1. / np.sqrt(N) * y0
        print('LNS (MEAN): (x0,y0) = ({},{}), (x1,y1) = ({},{})'.format(x0,y0,x1,y1))

        plt.legend(['Null', 'LNS (Max)', 'LNS (Mean)', 'Degenerate (Max)', 'Degenerate (Mean)', 'Ideal'])
        plt.ylabel('Std. Dev. [barns]', fontsize=16)
        plt.xlabel('# Histories [Millions]', fontsize=16)
        plt.xticks([1, 10, 100, 1000])
        plt.xlim([x0, x1])
        plt.savefig(subdirectory + filename, bbox_inches='tight')
