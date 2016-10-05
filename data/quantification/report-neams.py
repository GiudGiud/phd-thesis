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
os.chdir('reflector')

clusterizer_types = ['null', 'lns', 'pinch-agglomerative-(8)', 'degenerate']

min_capt = +1.e3
max_capt = -1.e3

for clusterizer_type in clusterizer_types:
    f = h5py.File('70-groups-{}.h5'.format(clusterizer_type))
    bias = f['70-groups'][clusterizer_type]['capture']['openmoc rel. err.']
    max_capt = max(np.nanmax(np.ravel(bias)), max_capt)
    min_capt = min(np.nanmin(np.ravel(bias)), min_capt)
    min_bias = np.nanmin(np.ravel(bias))
    max_bias = np.nanmax(np.ravel(bias))
    bias = max_bias if abs(max_bias) > abs(min_bias) else min_bias
    f.close()

if not os.path.exists('plots'):
    os.makedirs('plots')

# Instantiate Matplotlib figure with subplots
fig = plt.figure()
fig.set_figheight(8.5)
fig.set_figwidth(8.5)

cmap = plt.get_cmap('jet')
cmap.set_bad(alpha=0.0)

# Plot U-238 capture rate error heat maps
f = h5py.File('70-groups-null.h5')
bias = f['70-groups']['null']['capture']['openmoc rel. err.'][-1, ...]
f.close()

# plot a colormap of the fission rate percent rel. err.
plt.subplot(221)
im = plt.imshow(bias, interpolation='none',
                cmap=cmap, vmin=min_capt, vmax=max_capt)        
plt.grid(False)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlabel('Null', fontsize=16)

# Plot U-238 capture rate error heat maps
f = h5py.File('70-groups-lns.h5')
bias = f['70-groups']['lns']['capture']['openmoc rel. err.'][-1, ...]
f.close()

# plot a colormap of the fission rate percent rel. err.
plt.subplot(222)
im = plt.imshow(bias, interpolation='none',
                cmap=cmap, vmin=min_capt, vmax=max_capt)
plt.grid(False)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlabel('LNS', fontsize=16)

# Plot U-238 capture rate error heat maps
f = h5py.File('70-groups-pinch-agglomerative-(8).h5')
bias = f['70-groups']['pinch-agglomerative-(8)']['capture']['openmoc rel. err.'][-1, ...]
f.close()

# plot a colormap of the fission rate percent rel. err.
plt.subplot(223)
im = plt.imshow(bias, interpolation='none',
                cmap=cmap, vmin=min_capt, vmax=max_capt)
plt.grid(False)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlabel('Agglomerative', fontsize=16)

# Plot U-238 capture rate error heat maps
f = h5py.File('70-groups-degenerate.h5')
bias = f['70-groups']['degenerate']['capture']['openmoc rel. err.'][-1, ...]
f.close()

# plot a colormap of the fission rate percent rel. err.
plt.subplot(224)
im = plt.imshow(bias, interpolation='none',
                cmap=cmap, vmin=min_capt, vmax=max_capt)
plt.grid(False)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlabel('Degenerate', fontsize=16)

# Customize / shape the capture rates subplots
fig.subplots_adjust(bottom=0.15, hspace=0.15, wspace=0)
cbar_ax = fig.add_axes([0.165, 0.07, 0.7, 0.02])
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
for t in cbar.ax.get_xticklabels():
    t.set_fontsize(12)

plt.savefig('plots/capt-err.png', bbox_inches='tight', pad_inches=0)
