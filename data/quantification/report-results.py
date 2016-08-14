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
os.chdir(benchmark)

directories = {'assm-1.6': '1.6\\% Assm',
               'assm-3.1': '3.1\\% Assm',
               'assm-3.1-20BPs': '3.1\\% Assm w/ 20 BPs',
               '2x2': '2\\times2 Colorset',
               'reflector': '2\\times2 Colorset w/ Reflector',
               'full-core': 'BEAVRS Full Core'}

groups = [2, 2, 2] # 8, 70]
clusterizer_types = ['infinite', 'null', 'null'] #, 'degenerate']

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

print('MAX FISSION RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']['openmoc rel. err.']
        bias = np.nanmax(np.ravel(bias))
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
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']['openmoc rel. err.']
        bias = np.nanmean(np.ravel(bias))
        msg += '& {:1.2E} '.format(bias)
        f.close()
    msg += '\\\\\n'
print(msg)

print('MAX CAPTURE RATE ERROR')
msg = '\multirow{3}{*}{\parbox{2.5cm}{%s}} ' % directories[benchmark]
for clusterizer_type in clusterizer_types:
    msg += '& {} '.format(clusterizer_type.capitalize())
    for num_groups in groups:
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture']['openmoc rel. err.']
        bias = np.nanmax(np.ravel(bias))
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
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture']['openmoc rel. err.']
        bias = np.nanmean(np.ravel(bias))
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
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['fission']

        # Plot a colormap of the fission rate percent rel. err.
        subplot_ctr = '{}{}{}'.format(len(clusterizer_types),len(groups),i*len(groups)+j+1)
        plt.subplot(int(subplot_ctr))
        im = plt.imshow(
            bias['openmoc rel. err.'][-1, ...], interpolation='none', cmap=cmap)
        title = '{} ({}-group)'.format(clusterizer_type.capitalize(), num_groups)
        plt.title(title, fontsize=16)
        plt.grid(False)
        plt.axis('off')

        f.close()

# Customize / shape the fission rates subplots
fig.subplots_adjust(bottom=0.1, hspace=0.15, wspace=0)
cbar_ax = fig.add_axes([0.165, 0.07, 0.7, 0.02])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
for t in cbar.ax.get_xticklabels():
    t.set_fontsize(16)

plt.savefig('plots/fiss-err.png', bbox_inches='tight', pad_inches=0)

# Instantiate Matplotlib figure with subplots
fig = plt.figure()
fig.set_figheight(8.5)
fig.set_figwidth(8.5)

# Plot U-238 capture rate error heat maps
for i, clusterizer_type in enumerate(clusterizer_types):
    for j, num_groups in enumerate(groups):
        f = h5py.File('{}-groups-{}.h5'.format(num_groups, clusterizer_type))
        bias = f['{}-groups'.format(num_groups)][clusterizer_type]['capture']

        # Plot a colormap of the fission rate percent rel. err.
        subplot_ctr = '{}{}{}'.format(len(clusterizer_types),len(groups),i*len(groups)+j+1)
        plt.subplot(int(subplot_ctr))
        im = plt.imshow(
            bias['openmoc rel. err.'][-1, ...], interpolation='none', cmap=cmap)
        title = '{} ({}-groups)'.format(clusterizer_type.capitalize(), num_groups)
        plt.title(title, fontsize=16)
        plt.grid(False)
        plt.axis('off')

        f.close()

# Customize / shape the capture rates subplots
fig.subplots_adjust(bottom=0.1, hspace=0.15, wspace=0)
cbar_ax = fig.add_axes([0.165, 0.07, 0.7, 0.02])
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
for t in cbar.ax.get_xticklabels():
    t.set_fontsize(16)

plt.savefig('plots/capt-err.png', bbox_inches='tight', pad_inches=0)