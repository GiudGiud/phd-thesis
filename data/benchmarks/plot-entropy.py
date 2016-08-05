import os
import openmc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from _collections import OrderedDict


###############################################################################
# SINGLE FUEL ASSEMBLIES
###############################################################################

directories = OrderedDict()
directories['fuel-1.6'] = '1.6% Enr. (no BPs)'
directories['fuel-3.1'] = '3.1% Enr. (no BPs)'
directories['fuel-3.1-20BAs'] = '3.1% Enr. (20 BPs)'

# Create a matplotlib figure for all entropy convergence curves
fig = plt.figure()

for directory in directories:
    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.1000.h5'))
    batches = np.linspace(1, sp.current_batch, sp.current_batch, dtype=np.int)
    plt.plot(batches, sp.entropy, linewidth=2)

# Customize and save plot
plt.xlabel('Batch')
plt.ylabel('Shannon Entropy')
plt.legend(list(directories.values()), loc='center right')
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.4f'))
plt.savefig('entropy-assms.png', bbox_inches='tight')
plt.close()


###############################################################################
# MULTIPLE FUEL ASSEMBLIES
###############################################################################

# Multiple assembly directories
directories = OrderedDict()
directories['2x2'] = '2x2'
directories['reflector'] = 'reflector'
directories['full-core'] = 'full core'

# Create a matplotlib figure for all entropy convergence curves
fig = plt.figure()

for directory in directories:
    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.1000.h5'))
    batches = np.linspace(1, sp.current_batch, sp.current_batch, dtype=np.int)
    plt.plot(batches, sp.entropy/sp.entropy[-1], linewidth=2)

# Customize and save plot
plt.xlabel('Batch')
plt.ylabel('Shannon Entropy')
plt.legend(list(directories.values()))
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.savefig('entropy-multi-assms.png', bbox_inches='tight')
plt.close()

###############################################################################
# MULTIPLE FUEL ASSEMBLIES
###############################################################################

# Multiple assembly directories
directories = OrderedDict()
directories['fuel-1.6'] = '1.6% Enr. (no BPs)'
directories['fuel-3.1'] = '3.1% Enr. (no BPs)'
directories['fuel-3.1-20BAs'] = '3.1% Enr. (20 BPs)'
directories['2x2'] = '2x2 Colorset'
directories['reflector'] = '2x2 w/ Reflector'
directories['full-core'] = 'Full Core'

# Create a matplotlib figure for all entropy convergence curves
fig = plt.figure()

for directory in directories:
    sp = openmc.StatePoint(os.path.join(directory, 'statepoint.1000.h5'))
    batches = np.linspace(1, sp.current_batch, sp.current_batch, dtype=np.int)
    plt.plot(batches, sp.entropy/sp.entropy[-1], linewidth=2)

# Customize and save plot
plt.xlabel('Batch')
plt.ylabel('Normalized Shannon Entropy')
plt.legend(list(directories.values()))
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.savefig('entropy-all.png', bbox_inches='tight')
plt.close()
