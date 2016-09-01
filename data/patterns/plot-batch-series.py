import os

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


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict({'assm-1.6': '1.6% Enr. (no BPs)',
                           'assm-3.1':'3.1% Enr. (no BPs)',
                           'assm-3.1-20BPs': '3.1% Enr. (20 BPs)',
                           '2x2': '2x2 Colorset',
                           'reflector': '2x2 Colorset w/ Reflector',
                           'full-core': 'Full Core'})
directories = OrderedDict({'assm-1.6': '1.6% Enr. (no BPs)'})

coarse_groups = group_structures['CASMO']['2-group']

for directory in directories:
    print(directory)

    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.10000.h5'))

    # Load the MGXS library for this benchmark
    mgxslib = openmc.mgxs.Library.load_from_file(
        directory=os.path.join(directory, 'mgxs'), filename='distribcell')
    mgxslib.load_from_statepoint(sp)

    group = coarse_groups.get_group(6.67e-6)

    print('CAPTURE')

    # Loop over all fuel domains for each plot type
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))

        mgxs = mgxslib.get_mgxs(domain, 'capture')
        mgxs_pf = infermc.series.get_mgxs_batch_series(
            mgxs, coarse_groups=coarse_groups, xs_type='micro')

        print(mgxs_pf.shape)
        print(mgxs_pf[0].head(5))

        # Extract indices for the nuclide and energy of interest
        nuclide_indices = mgxs_pf.iloc[0]['nuclide'] == 'U238'
        energy_indices = mgxs_pf.iloc[0]['group'] == group

        # Slice the panel with the nuclide and energy indices
        pf_slice = mgxs_pf[:,nuclide_indices & energy_indices,:]

        # FIXME: Plot the max/mean absolute magnitude of the percent batch-by-batch change
        # FIXME: for each distribcell; as well as null; as well as LNS

        # FIXME: Make plot
