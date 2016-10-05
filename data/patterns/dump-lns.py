import os
import glob

import numpy as np
import pandas as pd
from _collections import OrderedDict

import openmc
import openmc.mgxs
import infermc
from infermc.energy_groups import group_structures


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict({'assm-1.6': '1.6% Enr. (no BPs)',
                           'assm-3.1':'3.1% Enr. (no BPs)',
                           'assm-3.1-20BPs': '3.1% Enr. (20 BPs)',
                           '2x2': '2x2 Colorset',
                           'reflector': '2x2 Colorset w/ Reflector'})
directories = OrderedDict({'full-core': 'Full Core'})

nuclide = 'U235'
mgxs_type = 'capture'
coarse_groups = group_structures['CASMO']['2-group']
cwd = os.getcwd()

for directory in directories:
    print(directory)

    os.chdir(os.path.join(cwd, directory))

    sp = openmc.StatePoint('statepoint.10000.h5')

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

        clusters = np.asarray(clusters)

        # Make filename
        if directory in ['2x2', 'reflector', 'full-core']:
            if '1.6' in domain.fill.name:
                filename = '16-enr.txt'
            elif '2.4' in domain.fill.name:
                filename = '24-enr.txt'
            elif '3.1' in domain.fill.name:
                filename = '31-enr.txt'
        else:
            filename = '{}.txt'.format(directory.replace('.', ''))

        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'lns')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass

        np.savetxt(subdirectory + filename, clusters)
