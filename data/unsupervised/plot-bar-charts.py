#!/usr/bin/env python

"""Plot Batchwise results including eigenvalue biases, reaction rates,
MGXS clusters, etc., by batch and clustering method."""


import os
import sys
import re
import glob
import shutil
import getopt
import random
from collections import OrderedDict
from distutils.util import strtobool

import h5py
import numpy as np
import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False                                                 # in your ~/.matplotlib/matplotlibrc

plt.ioff()

import matplotlib.pyplot as plt
import seaborn as sns
import infermc.plotter


# Initialize command line arguments
short_args = 'f:b:x:c:s:'
long_args = ['fix-cbar=', 'batches=', 'xscale=',
             'num-colors=', 'color-seed=', 'batchwise=']
opts, args = getopt.getopt(sys.argv[1:], short_args, long_args)

# Initialize default values for plotting parameters
fix_cbar = False
batches = 'last'
xscale = 'log'
num_colors = 8
color_seed = 20
batchwise = '70-groups' #'perfect'

# Parse command line arguments
for opt, arg in opts:
    if opt in ('-f', '--fix-cbar'):
        fix_cbar = strtobool(arg)
    elif opt in ('-b', '--batches'):
        batches = str(arg)
    elif opt in ('-x', '--xscale'):
        xscale = str(arg)
    elif opt in ('-c', '--num-colors'):
        num_colors = int(arg)
    elif opt in ('-s', '--color-seed'):
        color_seed = int(arg)
    elif opt in ('--batchwise'):
        batchwise = str(arg)

# Enter the directory of the benchmark of interest
benchmark = str(input('Benchmark: '))
os.chdir(benchmark)

# Create a directory for plots if it does not yet exist
if not os.path.exists('plots'):
    os.makedirs('plots')

# Open the file of Batchwise results
f = h5py.File('batchwise.h5', 'r')
batchwise = f[batchwise]

# Specify iterables of parameters over which to plot
clusterizers = ['openmc', 'null', 'degenerate', 'pinch-agglomerative-(4)',
                'pinch-agglomerative-(8)', 'pinch-agglomerative-(10)']
rxn_types = ['fission', 'capture']
metrics = ['max', 'mean']
nuclides = ['sum', 'U-238']

# Select randomized colors from seaborn
random.seed(color_seed)
colors = random.sample(list(sns.xkcd_rgb), num_colors)
cmap = sns.xkcd_palette(colors)

# Sort keys in "batchwise.h5" file such that legends are consistently ordered
batchwise_keys = \
    [key for key in batchwise if key not in ['null', 'degenerate', 'lns']]
batchwise_keys.insert(0, 'lns')
batchwise_keys.insert(0, 'degenerate')
batchwise_keys.insert(0, 'null')


################################################################################
#  Reaction Rate Error Bar Charts
################################################################################

# Identify batch index into HDF5 datasets for bar charts
if batches == 'first':
    batch = 0
else:
    batch = -1

# Plot the pin power relative error by the number of clusters
for rxn_type in rxn_types:
    for metric in metrics:
        for i, clusterizer in enumerate(clusterizers):

            # Extract relative error (uncertainty) for OpenMC
            num_clusters = ['openmc']
            rel_err = batchwise['null'][rxn_type]['openmc rel. err.'][batch, ...]
            if metric == 'max':
                data = [np.nanmax(np.fabs(rel_err))]
            elif metric == 'mean':
                data = [np.nanmean(np.fabs(rel_err))]

            num_clusters.append('degenerate')
            rel_err = batchwise['degenerate'][rxn_type]['openmoc rel. err.'][batch, ...]
            if metric == 'max':
                data.append(np.nanmax(np.fabs(rel_err)))
            elif metric == 'mean':
                data.append(np.nanmean(np.fabs(rel_err)))

            # Extract relative error for each clustering method
            for c in batchwise_keys:
                rel_err = batchwise[c][rxn_type]['openmoc rel. err.'][batch, ...]
                if metric == 'max':
                    rel_err = np.nanmax(np.abs(rel_err))
                else:
                    rel_err = np.nanmean(np.abs(rel_err))

                if clusterizer in c:
                    num_clusters.append(str(c.split('-')[-1][1:-1]))
                    data.append(rel_err)
                elif c == 'null':
                    num_clusters.append(c)
                    data.append(rel_err)

            # Plot a bar char for Clusterizer type
            ax = sns.barplot(x=num_clusters, y=data)
            ax.set_title('{}Clusterizer {} {} Rel. Err. [%]'.format(
                clusterizer.capitalize(), metric.capitalize(),
                rxn_type.capitalize()), fontsize=20)
            ax.set_xlabel('# clusters', fontsize=16)
            ax.set_ylabel('{} rel. err. [%]'.format(metric), fontsize=16)

            # Save the figure to a file or return to user if requested
            filename = 'plots/bar-{}-{}-rel-err-{}.png'.format(
                rxn_type, metric, clusterizer)
            plt.savefig(filename, bbox_inches='tight')
            plt.close(ax.figure)
