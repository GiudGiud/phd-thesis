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
clusterizers = ['pinch-agglomerative-(2)', 'pinch-agglomerative-(4)',
                'pinch-agglomerative-(6)', 'pinch-agglomerative-(8)']

# Select colors from seaborn
sns.palplot(sns.color_palette("hls", len(clusterizers)))

################################################################################
#  Plot Cluster Means
################################################################################

for clusterizer in clusterizers:
    for nuclide in batchwise[clusterizer]['clusters']['means']:
        for rxn_type in batchwise[clusterizer]['clusters']['means'][nuclide]:
            for group in batchwise[clusterizer]['clusters']['means'][nuclide][rxn_type]:
                for cell in batchwise[clusterizer]['clusters']['means'][nuclide][rxn_type][group]:
                    print(clusterizer, nuclide, group, cell)
                    means = batchwise[clusterizer]['clusters']['means'][nuclide][rxn_type][group][cell][...][:,:,0]
                    batches = batchwise[clusterizer]['batches'][...]
                
                    fig = plt.figure()
                    plt.semilogx(batches * num_particles, means, '-o')

                    # Annotate plot
                    title = '{} {} MGXS Cluster Means (Group {}/2)'.format(
                            nuclide, rxn_type.capitalize(), group.replace('group-', ''))
                    plt.title(title, fontsize=20)
                    plt.xlabel('# Histories', fontsize=16)
                    plt.ylabel('MGXS [barns]', fontsize=16)
                    plt.grid(True)
                    plt.ylim((0.00116, 0.00118))
                    
                    # Save the plot
                    filename = 'plots/evo-mgxs-mean-{}-{}-{}-{}.png'.format(
                        clusterizer, nuclide, rxn_type, group.replace('group-', ''))
                    plt.savefig(filename, bbox_inches='tight')
                    plt.close(fig)

################################################################################
#  Plot Cluster Population Variances
################################################################################

for clusterizer in clusterizers:
    for nuclide in batchwise[clusterizer]['clusters']['variances']:
        for rxn_type in batchwise[clusterizer]['clusters']['variances'][nuclide]:
            for group in batchwise[clusterizer]['clusters']['variances'][nuclide][rxn_type]:
                for cell in batchwise[clusterizer]['clusters']['variances'][nuclide][rxn_type][group]:
                    print(clusterizer, nuclide, group, cell)
                    variances = batchwise[clusterizer]['clusters']['variances'][nuclide][rxn_type][group][cell][...][:,:,0]
                    std_dev = np.sqrt(variances)
                    batches = batchwise[clusterizer]['batches'][...]
                    
                    fig = plt.figure()
                    plt.semilogx(batches * num_particles, std_dev, '-o')

                    # Annotate plot
                    title = '{} {} MGXS Cluster Std. Dev. (Group {}/2)'.format(
                        nuclide, rxn_type.capitalize(), group.replace('group-', ''))
                    plt.title(title, fontsize=20)
                    plt.xlabel('# Histories', fontsize=16)
                    plt.ylabel('Std. Dev. [barns]', fontsize=16)
                    plt.grid(True)
                    
                    # Save the plot
                    filename = 'plots/evo-mgxs-std-dev-{}-{}-{}-{}.png'.format(
                        clusterizer, nuclide, rxn_type, group.replace('group-', ''))
                    plt.savefig(filename, bbox_inches='tight')
                    plt.close(fig)

f.close()
