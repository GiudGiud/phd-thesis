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
quarter_symmetric = str(input('Quarter Core Symmetric: '))
code = str(input('openmc or openmoc: '))
metric = str(input('mean or rel. err.: '))

if metric == 'rel. err.':
    hdf5_key = '{} {}'.format(code.lower(), metric.lower())
else:
    hdf5_key = code.lower()
if quarter_symmetric == 'True':
    quarter_symmetric = True
else:
    quarter_symmetric = False

os.chdir(benchmark)

directories = {'assm-1.6': '1.6\\% Assm',
               'assm-3.1': '3.1\\% Assm',
               'assm-3.1-20BPs': '3.1\\% Assm w/ 20 BPs',
               '2x2': '2$\\times$2 Colorset',
               'reflector': '2$\\times$2 Colorset w/ Reflector',
               'full-core': 'BEAVRS Full Core'}

groups = [2] #, 8, 70]
clusterizer_types = ['infinite', 'null', 'degenerate']

print('EIGENVALUE BIAS')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['keff']['bias']
        msg += '& {:.0f} '.format(bias[0])
        f.close()
    msg += '\\\\\n'
print(msg)

min_fiss = +1.e3
max_fiss = -1.e3

print('MAX FISSION RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission'][hdf5_key]
        max_fiss = max(np.nanmax(np.ravel(bias)), max_fiss)
        min_fiss = min(np.nanmin(np.ravel(bias)), min_fiss)
        min_bias = np.nanmin(np.ravel(bias))
        max_bias = np.nanmax(np.ravel(bias))
        bias = max_bias if abs(max_bias) > abs(min_bias) else min_bias
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
print(msg)

print('MEAN FISSION RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission'][hdf5_key]
        bias = np.nanmean(np.fabs(np.ravel(bias)))
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
print(msg)

min_capt = +1.e3
max_capt = -1.e3

print('MAX CAPTURE RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture'][hdf5_key]
        max_capt = max(np.nanmax(np.ravel(bias)), max_fiss)
        min_capt = min(np.nanmin(np.ravel(bias)), min_fiss)
        min_bias = np.nanmin(np.ravel(bias))
        max_bias = np.nanmax(np.ravel(bias))
        bias = max_bias if abs(max_bias) > abs(min_bias) else min_bias
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
print(msg)

print('MEAN CAPTURE RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture'][hdf5_key]
        bias = np.nanmean(np.fabs(np.ravel(bias)))
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
print(msg)

if not os.path.exists('plots'):
    os.makedirs('plots')

# Instantiate Matplotlib figure with subplots
fig = plt.figure()
fig.set_figheight(8.5)
fig.set_figwidth(8.5)

cmap = plt.get_cmap('jet')
cmap.set_bad(alpha=0.0)

# Plot fission rate error heat maps
for i, clusterizer_type in enumerate(clusterizer_types):
    for j, num_groups in enumerate(groups):
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission'][hdf5_key][-1, ...]

        # Extract the array and make it quarter core symmetric
        if quarter_symmetric:
            full_bias = np.zeros((2*bias.shape[0]-1, 2*bias.shape[0]-1), dtype=np.float)
            full_bias[:bias.shape[0], :bias.shape[1]] = bias[::-1,:]    # top left
            full_bias[bias.shape[0]-1:, :bias.shape[1]] = bias[:,:]       # top right
            full_bias[:bias.shape[0], bias.shape[1]-1:] = bias[::-1,::-1]     # bottom left
            full_bias[bias.shape[0]-1:, bias.shape[1]-1:] = bias[:,::-1]     # bottom right

            if metric == 'mean':
                full_bias[bias.shape[0]-1,:] *= 2.
                full_bias[:,bias.shape[1]-1] *= 2.

            bias = full_bias

        # Plot a colormap of the fission rate percent rel. err.
        subplot_ctr = '{}{}{}'.format(len(clusterizer_types),len(groups),i*len(groups)+j+1)
        plt.subplot(int(subplot_ctr))
        im = plt.imshow(bias, interpolation='none',
            cmap=cmap, vmin=min_fiss, vmax=max_fiss)

        if num_groups == groups[0]:
            plt.ylabel(clusterizer_type.capitalize(), fontsize=16)
        if clusterizer_type == clusterizer_types[0]:
            plt.title('{} Groups'.format(num_groups), fontsize=16)
        plt.grid(False)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)

        f.close()

# Customize / shape the fission rates subplots
fig.subplots_adjust(bottom=0.1, hspace=0.15, wspace=0)
cbar_ax = fig.add_axes([0.165, 0.07, 0.7, 0.02])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
for t in cbar.ax.get_xticklabels():
    t.set_fontsize(12)

plt.savefig('plots/fiss-err.png', bbox_inches='tight', pad_inches=0)

# Instantiate Matplotlib figure with subplots
fig = plt.figure()
fig.set_figheight(8.5)
fig.set_figwidth(8.5)

# Plot U-238 capture rate error heat maps
for i, clusterizer_type in enumerate(clusterizer_types):
    for j, num_groups in enumerate(groups):
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture'][hdf5_key][-1, ...]

        # Extract the array and make it quarter core symmetric
        if quarter_symmetric:
            full_bias = np.zeros((2*bias.shape[0]-1, 2*bias.shape[0]-1), dtype=np.float)
            full_bias[:bias.shape[0], :bias.shape[1]] = bias[::-1,:]    # top left
            full_bias[bias.shape[0]-1:, :bias.shape[1]] = bias[:,:]       # top right
            full_bias[:bias.shape[0], bias.shape[1]-1:] = bias[::-1,::-1]     # bottom left
            full_bias[bias.shape[0]-1:, bias.shape[1]-1:] = bias[:,::-1]     # bottom right

            if metric == 'mean':
                full_bias[bias.shape[0]-1,:] *= 2.
                full_bias[:,bias.shape[1]-1] *= 2.

            bias = full_bias

        # Plot a colormap of the fission rate percent rel. err.
        subplot_ctr = '{}{}{}'.format(len(clusterizer_types),len(groups),i*len(groups)+j+1)
        plt.subplot(int(subplot_ctr))
        im = plt.imshow(bias, interpolation='none',
            cmap=cmap, vmin=min_fiss, vmax=max_fiss)

        if num_groups == groups[0]:
            plt.ylabel(clusterizer_type.capitalize(), fontsize=16)
        if clusterizer_type == clusterizer_types[0]:
            plt.title('{} Groups'.format(num_groups), fontsize=16)
        plt.grid(False)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)

        f.close()

# Customize / shape the capture rates subplots
fig.subplots_adjust(bottom=0.1, hspace=0.15, wspace=0)
cbar_ax = fig.add_axes([0.165, 0.07, 0.7, 0.02])
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
for t in cbar.ax.get_xticklabels():
    t.set_fontsize(12)

plt.savefig('plots/capt-err.png', bbox_inches='tight', pad_inches=0)
