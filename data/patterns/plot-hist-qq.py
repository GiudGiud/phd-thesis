import os

import matplotlib

# force headless backend, or set 'backend' to 'Agg'                                                          # in your ~/.matplotlib/matplotlibrc                                                                         
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
import sklearn.preprocessing
from _collections import OrderedDict

import openmc
from infermc.energy_groups import group_structures
from infermc.plotter import hist_kde_rug, normal_qq


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict({'assm-1.6': '1.6% Enr. (no BPs)',
                           'assm-3.1':'3.1% Enr. (no BPs)',
                           'assm-3.1-20BPs': '3.1% Enr. (20 BPs)',
                           '2x2': '2x2 Colorset',
                           'reflector': '2x2 Colorset w/ Reflector',
                           'full-core': 'Full Core'})

for directory in directories:
    print(directory)
    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.10000.h5'))
    batches = np.linspace(1, sp.current_batch, sp.current_batch, dtype=np.int)

    # Load the MGXS library for this benchmark
    coarse_groups = group_structures['CASMO']['2-group']
    mgxslib = openmc.mgxs.Library.load_from_file(
        directory=os.path.join(directory, 'mgxs'), filename='distribcell')
    mgxslib.load_from_statepoint(sp)
    mgxslib = mgxslib.get_condensed_library(coarse_groups)

    group = coarse_groups.get_group(6.67e-6)

    print('CAPTURE')

    # FIXME: Loop over all fuel domains for each plot type
    # FIXME: Plot histogram-kde-rug plot
    # FIXME: QQ plot
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))
        mgxs = mgxslib.get_mgxs(domain, 'capture')
        mgxs_df = mgxs.get_pandas_dataframe(
            xs_type='micro', nuclides=['U238'], 
            groups=[1], distribcell_paths=False)
        mgxs_df = mgxs_df[mgxs_df['mean'] > 0]
        hist_kde_rug(mgxs_df, domain, 'U238', 'capture', 1, directory)
        normal_qq(mgxs_df, domain, 'U238', 'capture', 1, False, directory)

    print('FISSION')

    # FIXME: Loop over all fuel domains for each plot type
    # FIXME: Plot histogram-kde-rug plot
    # FIXME: QQ plot             
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))
        mgxs = mgxslib.get_mgxs(domain, 'fission')
        mgxs_df = mgxs.get_pandas_dataframe(
            xs_type='micro', nuclides=['U235'],
            groups=[2], distribcell_paths=False)
        mgxs_df = mgxs_df[mgxs_df['mean'] > 0]
        hist_kde_rug(mgxs_df, domain, 'U235', 'fission', 2, directory)
        normal_qq(mgxs_df, domain, 'U235', 'fission', 2, False, directory)
