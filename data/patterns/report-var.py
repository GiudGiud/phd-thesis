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

directories = OrderedDict()
directories['assm-1.6-inf'] = '1.6% Enr. (no \acp{CRGT}, no \acp{BP})'
directories['assm-3.1-inf'] = '3.1% Enr. (no \acp{CRGT}, no \acp{BP})'
directories['assm-1.6'] = '1.6% Enr. (no BPs)'
directories['assm-3.1'] = '3.1% Enr. (no BPs)'
directories['assm-3.1-20BPs'] = '3.1% Enr. (20 BPs)'
directories['2x2'] = '2x2 Colorset'
directories['reflector'] = '2x2 Colorset w/ Reflector'
directories['full-core'] = 'Full Core'

print('CAPTURE')

table16 = ''
table31 = ''

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

    flag16 = False
    flag31 = False

    # Loop over all fuel domains for each plot type
    for domain in mgxslib.domains:
        mgxs = mgxslib.get_mgxs(domain, 'capture')
        mgxs_df = mgxs.get_pandas_dataframe(
            xs_type='micro', nuclides=['U238'], 
            groups=[1], distribcell_paths=False)
        mgxs_df = mgxs_df[mgxs_df['std. dev.'] > 0]

        if '1.6' in domain.fill.name:
            table16 += directories[directory]
            print(mgxs_df['mean'].abs().max(), mgxs_df['mean'].abs().mean())
            table16 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().max())
            table16 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().mean())
            flag16 = True
        elif '3.1' in domain.fill.name:
            table31 += directories[directory]
            table31 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().max())
            table31 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().mean())
            flag31 = True

    if flag16:
        table16 += '\midrule\n'
        flag16 = False
    if flag31:
        table31 += '\midrule\n'
        flag31 = False

table = '\multicolumn{3}{c}{\bf 1.6\% Enr. Fuel} \\\\\n'
table += '\midrule\n' + table16 + '\midrule\n'
table += '\multicolumn{3}{c}{\bf 3.1\% Enr. Fuel} \\\\\n'
table += table31 + '\\bottomrule\n'
print(table)
print()

print('FISSION')

table16 = ''
table31 = ''

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

    group = coarse_groups.get_group(6.67e-7)

    flag16 = False
    flag31 = False

    # Loop over all fuel domains for each plot type
    for domain in mgxslib.domains:
        mgxs = mgxslib.get_mgxs(domain, 'total')
        mgxs_df = mgxs.xs_tally.get_pandas_dataframe()
        mgxs_df = mgxs.get_pandas_dataframe(
            xs_type='micro', nuclides=['U235'],
            groups=[2], distribcell_paths=False)
        mgxs_df = mgxs_df[mgxs_df['std. dev.'] > 0]

        if '1.6' in domain.fill.name:
            table16 += directories[directory]
            table16 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().max())
            table16 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().mean())
            flag16 = True
        elif '3.1' in domain.fill.name:
            table31 += directories[directory]
            table31 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().max())
            table31 += ' & {:1.3E} \\\\\n'.format(mgxs_df['std. dev.'].abs().mean())
            flag31 = True

    if flag16:
        table16 += '\midrule\n'
        flag16 = False
    if flag31:
        table31 += '\midrule\n'
        flag31 = False

table = '\multicolumn{3}{c}{\bf 1.6\% Enr. Fuel} \\\\\n'
table += '\midrule\n' + table16 + '\midrule\n'
table += '\multicolumn{3}{c}{\bf 3.1\% Enr. Fuel} \\\\\n'
table += table31 + '\\bottomrule\n'
print(table)
print()
