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

    # Loop over all fuel domains for each plot type
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))
        mgxs = mgxslib.get_mgxs(domain, 'capture')
        mgxs_df = mgxs.get_pandas_dataframe(
            xs_type='micro', nuclides=['U238'], 
            groups=[1], distribcell_paths=False)
        mgxs_df = mgxs_df[mgxs_df['mean'] > 0]

        # Make title
        suptitle = 'U-238 Capture MGXS (Group 1/2)'
        title = directories[directory]
        if directory in ['2x2', 'reflector', 'full-core']:
            if '1.6' in domain.fill.name:
                title += ' - 1.6% Enr. Fuel'
                filename = '16-enr-capt-1.png'
            elif '2.4' in domain.fill.name:
                title += ' - 2.4% Enr. Fuel'
                filename = '24-enr-capt-1.png'
            elif '3.1' in domain.fill.name:
                title += ' - 3.1% Enr. Fuel'
                filename = '31-enr-capt-1.png'
        else:
            filename = '{}-capt-1.png'.format(directory.replace('.', ''))

        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'hist-kde-rug/')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass

        # Make filename
#        filename = subdirectory + '{}-capt-1.png'.format(directory)

        fig = hist_kde_rug(mgxs_df, domain, 'U238',
                           'capture', 1, directory, get_figure=True)
        axes = fig.get_axes()[0]
        axes.set_title('', fontsize=16)
        fig.suptitle(title, fontsize=20)
        axes.tick_params(axis='both', which='major', labelsize=16)
        axes.tick_params(axis='both', which='minor', labelsize=16)
        axes.set_xlabel('Mean [barns]', fontsize=16)
        axes.set_ylabel('# Samples', fontsize=16)
        fig.savefig(subdirectory + filename, bbox_inches='tight')
        plt.close(fig)

        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'quantile/')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass

        # Make filename
#        filename = subdirectory + '{}-qq-capt-1.png'.format(directory)

        fig = normal_qq(mgxs_df, domain, 'U238',
                        'capture', 1, True, directory, get_figure=True)
        fig = fig.figure
        axes = fig.get_axes()[0]
        axes.set_title('', fontsize=20)
        fig.suptitle(title, fontsize=20)
        axes.tick_params(axis='both', which='major', labelsize=16)
        axes.tick_params(axis='both', which='minor', labelsize=16)
        axes.set_xlabel('Theoretical Quantiles', fontsize=16)
        axes.set_ylabel('Sample Quantiles', fontsize=16)
        fig.savefig(subdirectory + filename, bbox_inches='tight')
        plt.close(fig)

    print('FISSION')

    # Loop over all fuel domains for each plot type           
    for domain in mgxslib.domains:
        print('cell {}'.format(domain.id))
        mgxs = mgxslib.get_mgxs(domain, 'fission')
        mgxs_df = mgxs.get_pandas_dataframe(
            xs_type='micro', nuclides=['U235'],
            groups=[2], distribcell_paths=False)
        mgxs_df = mgxs_df[mgxs_df['mean'] > 0]

        # Make title
        suptitle = 'U-235 Fission MGXS (Group 2/2)'
        title = directories[directory]
        if directory in ['2x2', 'reflector', 'full-core']:
            if '1.6' in domain.fill.name:
                title += ' - 1.6% Enr. Fuel'
                filename = '16-enr-fiss-2.png'
            elif '2.4' in domain.fill.name:
                title += ' - 2.4% Enr. Fuel'
                filename = '24-enr-fiss-2.png'
            elif '3.1' in domain.fill.name:
                title += ' - 3.1% Enr. Fuel'
                filename = '31-enr-fiss-2.png'
        else:
            filename = '{}-fiss-2.png'.format(directory.replace('.', ''))

        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'hist-kde-rug/')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass

        # Make filename
#        filename = subdirectory + '{}-fiss-2.png'.format(directory)

        fig = hist_kde_rug(mgxs_df, domain, 'U235',
                           'fission', 2, directory, get_figure=True)
        axes = fig.get_axes()[0]
        axes.set_title('', fontsize=16)
        fig.suptitle(title, fontsize=20)
        axes.tick_params(axis='both', which='major', labelsize=16)
        axes.tick_params(axis='both', which='minor', labelsize=16)
        axes.set_xlabel('Mean [barns]', fontsize=16)
        axes.set_ylabel('# Samples', fontsize=16)
        fig.savefig(subdirectory + filename, bbox_inches='tight')
        plt.close(fig)

        # Make directory if it does not exist
        subdirectory = os.path.join('plots', directory, 'quantile/')
        try:
            os.makedirs(subdirectory)
        except OSError:
            pass

        # Make filename
#        filename = subdirectory + '{}-fiss-2.png'.format(directory)

        fig = normal_qq(mgxs_df, domain, 'U235',
                        'fission', 2, True, directory, get_figure=True)
        fig = fig.figure
        axes = fig.get_axes()[0]
        axes.set_title('', fontsize=16)
        fig.suptitle(title, fontsize=20)
        axes.tick_params(axis='both', which='major', labelsize=16)
        axes.tick_params(axis='both', which='minor', labelsize=16)
        axes.set_xlabel('Theoretical Quantiles', fontsize=16)
        axes.set_ylabel('Sample Quantiles', fontsize=16)
        fig.savefig(subdirectory + filename, bbox_inches='tight')
        plt.close(fig)
