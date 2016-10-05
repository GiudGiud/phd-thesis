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
import seaborn as sns
import matplotlib.pyplot as plt

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
batchwise = 'perfect'

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

# Create a directory for plots if it does not yet exist
if not os.path.exists('plots'):
    os.makedirs('plots')

# Open the file of Batchwise results
f = h5py.File('batchwise.h5', 'r')
batchwise = f[batchwise]

# Specify iterables of parameters over which to plot
clusterizers = ['pinch', 'combined', 'local']
rxn_types = ['fission', 'capture']
metrics = ['max', 'mean']
nuclides = ['sum', 'U-238']

# Select randomized colors from seaborn
random.seed(color_seed)
colors = random.sample(list(sns.xkcd_rgb), num_colors)
cmap = sns.xkcd_palette(colors)

# Sort keys in "batchwise.h5" file such that legends are consistently ordered
batchwise_keys = \
    [key for key in batchwise if key not in ['null', 'degenerate', 'infinite']]
batchwise_keys.insert(0, 'degenerate')
batchwise_keys.insert(0, 'null')
batchwise_keys.insert(0, 'infinite')


################################################################################
#  Eigenvalue Bias by Batch
################################################################################

# Plot keff bias separately for pinch, combined and local clustering
for clusterizer in clusterizers:
    pattern = r'null|degenerate|infinite|{}'.format(clusterizer)
    regex = re.compile(pattern)
    fig = plt.figure()
    legend = []
    color = 0

    # Add each clustering
    for c in batchwise_keys:
        if regex.match(c) is not None:
            if xscale == 'linear':
                plt.plot(batchwise[c]['batches'][...],
                         batchwise[c]['keff']['bias'],
                         '-o', color=cmap[color])
            elif xscale == 'log':
                plt.semilogx(batchwise[c]['batches'][...],
                             batchwise[c]['keff']['bias'],
                             '-o', color=cmap[color])
            legend.append(c)
            color += 1

    # Annotate plot
    plt.title('Eigenvalue Bias', fontsize=20)
    plt.xlabel('batch', fontsize=16)
    plt.ylabel('bias [pcm]', fontsize=16)
    plt.legend(legend, loc='best')
    plt.grid(True)
    filename = 'plots/keff-bias-{}.png'.format(clusterizer)
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


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

            num_clusters.append('infinite')
            rel_err = batchwise['infinite'][rxn_type]['openmoc rel. err.'][batch, ...]
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


################################################################################
#  Reaction Rate Errors by Batch
################################################################################

def rel_err_stats(batchwise, rxn_type='fission',
                  metric='max', exclude=None):
    """Compute the relative error of various OpenMOC reaction rates
    with respect to fully converged OpenMC results.

    Parameters
    ----------
    batchwise : h5py.Group
        An HDF5 group of Batchwise results with reaction rates
    rxn_type : {'fission', 'capture'}
        Reaction rate of interest
    metric : {'max', 'mean'}
        The metric to compute for each batch
    exclude : re.SRE_Pattern or None, optional
        Exclude clustering algorithms which match this regular expression in
        their key in the HDF5 file

    Returns
    -------
    rel_err : dict
        A dictionary of the batchwise reaction rate relative errors
        indexed by clustering algorithm.

    """

    # Initialize a container for the relative errors
    rel_err = OrderedDict()

    # Store the percent relative error for each clustering algorithm
    # which does not meet the exclude criteria
    for c in batchwise_keys:
        if exclude is None or exclude.match(c) is None:
            rel_err[c] = \
                np.fabs(batchwise[c][rxn_type]['openmoc rel. err.'][...])
            rel_err['openmc'] = \
                np.fabs(batchwise[c][rxn_type]['openmc rel. err.'][...])

            if metric == 'max':
                rel_err[c] = np.nanmax(rel_err[c], axis=(1,2))
                rel_err['openmc'] = np.nanmax(rel_err['openmc'], axis=(1,2))
            elif metric == 'mean':
                rel_err[c] = np.nanmean(rel_err[c], axis=(1,2))
                rel_err['openmc'] = np.nanmean(rel_err['openmc'], axis=(1,2))

    return rel_err


# Plot the pin power relative error evolution
for i, clusterizer in enumerate(clusterizers):
    exclude = clusterizers[:i] + clusterizers[(i+1):]
    exclude = r'lns|{}|{}'.format(exclude[0], exclude[1])
    exclude = re.compile(exclude)

    for rxn_type, nuclide in zip(rxn_types, nuclides):
        for metric in metrics:
            rel_err = rel_err_stats(batchwise, rxn_type, metric, exclude)
            fig = plt.figure()
            color = 0

            # Plot the relative errors for each method under
            # this family of clustering algorithms
            for c in rel_err:
                if xscale == 'linear':
                    plt.plot(batchwise['null']['batches'][...],
                             rel_err[c], '-o', color=cmap[color])
                elif xscale == 'log':
                    plt.semilogx(batchwise['null']['batches'][...],
                                 rel_err[c], '-o', color=cmap[color])
                color += 1

            # Annotate plot
            title = '{} {} Rate Rel. Err.'.format(
                metric.capitalize(), rxn_type.capitalize())
            plt.title(title, fontsize=20)
            plt.xlabel('batch', fontsize=16)
            plt.ylabel('{} error [%]'.format(metric), fontsize=16)
            plt.legend(rel_err.keys(), loc='best')
            plt.grid(True)

            # Save the plot
            filename = 'plots/evo-{}-{}-{}.png'.format(rxn_type, metric, clusterizer)
            plt.savefig(filename, bbox_inches='tight')
            plt.close(fig)


################################################################################
#  Reaction Rate Error Distribution Histograms
################################################################################

# Plot a histogram of the reaction rate relative errors
for c in batchwise_keys:
    for rxn_type in rxn_types:
        rxn_rates = batchwise[c][rxn_type]['openmoc rel. err.']
        rxn_rates = rxn_rates[-1, ...].flatten()
        rxn_rates = rxn_rates[~np.isnan(rxn_rates)]
        ax = sns.distplot(rxn_rates, rug=True, norm_hist=False,
                          rug_kws={"color": "g"},
                          kde_kws={"color": "k", "lw": 3, "label": "KDE"},
                          hist_kws={"histtype": "step", "linewidth": 3,
                                  "alpha": 1, "color": "b"})
        ax.set_title('{} {} Rel. Err. [%]'.format(
            c.capitalize(), rxn_type.capitalize()))
        ax.set_xlabel('rel. err. [%]')
        ax.set_ylabel('# samples')

        # Save the figure to a file or return to user if requested
        filename = 'plots/hist-rel-err-{}-{}.png'.format(rxn_type, c)
        plt.savefig(filename, bbox_inches='tight')
        plt.close(ax.figure)
