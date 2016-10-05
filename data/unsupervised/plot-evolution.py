#!/usr/bin/env python

"""Plot Batchwise results including eigenvalue biases, reaction rates,
MGXS clusters, etc., by batch and clustering method."""


import os
import glob
from collections import OrderedDict

import h5py
import numpy as np
import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc

plt.ioff()

import matplotlib.pyplot as plt
import seaborn as sns
import infermc.plotter

# Initialize default values for plotting parameters
batchwise = '70-groups'
num_particles = 1E5 # particles / batch

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
clusterizers = ['openmc', 'null', 'degenerate', 'pinch-kmeans-(4)',
                'pinch-kmeans-(8)', 'pinch-kmeans-(10)']
#clusterizers = ['pinch']
rxn_types = ['fission', 'capture']
metrics = ['max', 'mean']
nuclides = ['sum', 'U-238']

# Select colors from seaborn
sns.palplot(sns.color_palette("hls", len(clusterizers)))

# Sort keys in "batchwise.h5" file such that legends are consistently ordered
batchwise_keys = \
    [key for key in batchwise if key not in ['null', 'degenerate', 'lns']]
batchwise_keys.insert(0, 'lns')
batchwise_keys.insert(0, 'degenerate')
batchwise_keys.insert(0, 'null')


################################################################################
#  Eigenvalue Bias by Batch
################################################################################

fig = plt.figure()
legend = []

# Plot keff bias separately for pinch, combined and local clustering
for c in clusterizers:
    if c == 'openmc':
        continue

    # Add each clustering
    plt.semilogx(batchwise[c]['batches'][...] * num_particles,
                 batchwise[c]['keff']['bias'], '-o')

    if 'pinch-agglomerative' in c:
        legend.append(c.replace('pinch-agglomerative', 'imgxs'))
    elif 'pinch-agglomerative' in c:
        legend.append(c.replace('combined-agglomerative', 'imgxs'))
    else:
        legend.append(c)

# Annotate plot
plt.title('Eigenvalue Bias', fontsize=20)
plt.xlabel('# Histories', fontsize=16)
plt.ylabel(r'Bias $\Delta\rho$ [pcm]', fontsize=16)
plt.legend(legend, loc='best')
plt.grid(True)
filename = 'plots/keff-bias.png'
plt.savefig(filename, bbox_inches='tight')
plt.close(fig)


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


for rxn_type, nuclide in zip(rxn_types, nuclides):
    for metric in metrics:

        rel_err = rel_err_stats(batchwise, rxn_type, metric)
        fig = plt.figure()
        legend = []

        for c in clusterizers:

            if c == 'openmc':
                start = 1
            else:
                start = 0
                
            plt.semilogx(
                batchwise['null']['batches'][...][start:] * num_particles, rel_err[c][start:], '-o')

            if 'pinch-agglomerative' in c:
                legend.append(c.replace('pinch-agglomerative', 'imgxs'))
            elif 'pinch-agglomerative' in c:
                legend.append(c.replace('combined-agglomerative', 'imgxs'))
            else:
                legend.append(c)

        # Annotate plot
        title = '{} {} Rate Rel. Err.'.format(
            metric.capitalize(), rxn_type.capitalize())
        plt.title(title, fontsize=20)
        plt.xlabel('# Histories', fontsize=16)
        plt.ylabel('Rel. Err. [%]'.format(metric), fontsize=16)
        plt.legend(legend, loc='best')
        plt.grid(True)
        
        # Save the plot
        filename = 'plots/evo-{}-{}.png'.format(rxn_type, metric)
        plt.savefig(filename, bbox_inches='tight')
        plt.close(fig)
