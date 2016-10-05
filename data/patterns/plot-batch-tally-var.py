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
        group = coarse_groups.get_group(0.2)
    else:
        group = coarse_groups.get_group(1.e-9)
    lower_bound = coarse_groups.get_group_bounds(group)[0]

    # Loop over all fuel domains for each plot type
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))

        mgxs = mgxslib.get_mgxs(domain, mgxs_type)
        tally = mgxs.tallies['flux']
        std_devs = []
        batches = []

        for filename in statepoints:
            print(filename)
            sp = openmc.StatePoint(filename)
            if sp.tallies_present:
                tally = sp.get_tally(id=tally.id)
                batches.append(sp.n_realizations)
                std_devs.append(tally.std_dev.ravel()[0])

        print(batches)
        print(std_devs)

        if '1.6' in domain.fill.name:
            filename = '16-enr-var-flux-{}.png'.format(group)
        elif '2.4' in domain.fill.name:
            filename = '24-enr-var-flux-{}.png'.format(group)
        elif '3.1' in domain.fill.name:
            filename = '31-enr-var-flux-{}.png'.format(group)
        else:
            filename = '{}-var-flux-{}.png'.format(directory.replace('.', ''), group)

        print('here')
        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'convergence/')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass


        print('creating figures')
        # Plot the figure
        plt.figure()
        plt.loglog(batches, std_devs, linestyle='-', color='b')

        # FIXME: Line for NULL
        N = 9900 - 20
        x0 = batches[0]
        x1 = batches[-1]
        y0 = std_devs[0]
        y1 = 1. / np.sqrt(N) * y0
        print('NULL: (x0,y0) = ({},{}), (x1,y1) = ({},{})'.format(x0,y0,x1,y1))
        plt.loglog([x0, x1], [y0, y1], linestyle='-', color='g', linewidth=0.5)

        plt.legend(['Observed', 'Ideal'])
        plt.ylabel('Std. Dev.', fontsize=16)
        plt.xlabel('# Batches', fontsize=16)
#        plt.xticks([1, 10, 100, 1000])
#        plt.xlim([1, 1000])
        plt.savefig(subdirectory + filename, bbox_inches='tight')
