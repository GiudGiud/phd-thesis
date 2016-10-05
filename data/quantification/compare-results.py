import os

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

# Enter the directory of the benchmark of interest
benchmark = str(input('Benchmark: '))
reference = str(input('Reference: '))
evaluate = str(input('Evaluate: '))
num_groups = 70

reference_f = h5py.File(os.path.join(benchmark, '{}-groups-{}.h5'.format(num_groups, reference)))
evaluate_f = h5py.File(os.path.join(benchmark, '{}-groups-{}.h5'.format(num_groups, evaluate)))

reference_rates = reference_f['{}-groups'.format(num_groups)][reference]['fission']['openmoc'][-1, ...]
evaluate_rates = evaluate_f['{}-groups'.format(num_groups)][evaluate]['fission']['openmoc'][-1, ...]
delta_rates = (reference_rates - evaluate_rates) / reference_rates * 100.
delta_rates = delta_rates[::-1, ::-1]

reference_f.close()
evaluate_f.close()

if not os.path.exists(os.path.join(benchmark, 'plots')):
    os.makedirs(os.path.join(benchmark, 'plots'))

cmap = plt.get_cmap('jet')
cmap.set_bad(alpha=0.0)

# Instantiate Matplotlib figure with subplots
fig = plt.figure()
im = plt.imshow(delta_rates, interpolation='none', cmap=cmap)
plt.title('{}-to-{} Fission % ({} Groups)'.format(reference, evaluate, num_groups), fontsize=16)
plt.colorbar()
plt.grid(False)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.savefig(os.path.join(benchmark, 'plots', 'fiss-comp.png'), bbox_inches='tight', pad_inches=0)



reference_f = h5py.File(os.path.join(benchmark, '{}-groups-{}.h5'.format(num_groups, reference)))
evaluate_f = h5py.File(os.path.join(benchmark, '{}-groups-{}.h5'.format(num_groups, evaluate)))

reference_rates = reference_f['{}-groups'.format(num_groups)][reference]['capture']['openmoc'][-1, ...]
evaluate_rates = evaluate_f['{}-groups'.format(num_groups)][evaluate]['capture']['openmoc'][-1, ...]
delta_rates = (reference_rates - evaluate_rates) / reference_rates * 100.
delta_rates = delta_rates[::-1, ::-1]

reference_f.close()
evaluate_f.close()

if not os.path.exists(os.path.join(benchmark, 'plots')):
    os.makedirs(os.path.join(benchmark, 'plots'))

cmap = plt.get_cmap('jet')
cmap.set_bad(alpha=0.0)

# Instantiate Matplotlib figure with subplots
fig = plt.figure()
im = plt.imshow(delta_rates, interpolation='none', cmap=cmap)
plt.title('{}-to-{} U-238 Capture % ({} Groups)'.format(reference, evaluate, num_groups), fontsize=16)
plt.colorbar()
plt.grid(False)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.savefig(os.path.join(benchmark, 'plots', 'capt-comp.png'), bbox_inches='tight', pad_inches=0)
