import openmc.mgxs
from infermc.energy_groups import group_structures

print('anisotropic')

# Load MGXS Library and condense to 1 group
mgxs_lib = openmc.mgxs.Library.load_from_file(directory='anisotropic')
coarse_groups = group_structures['CASMO']['{}-group'.format(1)]
mgxs_lib = mgxs_lib.get_condensed_library(coarse_groups)

# Extract MGXS objects from library
nufiss = mgxs_lib.get_mgxs(mgxs_lib.domains[0], 'nu-fission')
total = mgxs_lib.get_mgxs(mgxs_lib.domains[0], 'total')
scatt = mgxs_lib.get_mgxs(mgxs_lib.domains[0], 'nu-scatter matrix')

# Get NumPy arrays of data
nufiss = nufiss.get_xs(nuclides='sum')
total = total.get_xs(nuclides='sum')
scatt = scatt.get_xs(nuclides='sum')
absorb = total - scatt

keff = nufiss / absorb
print('keff = {0:1.6f}'.format(float(keff[0])))


print('iso-in-lab')

# Load MGXS Library and condense to 1 group
mgxs_lib = openmc.mgxs.Library.load_from_file(directory='iso-in-lab')
coarse_groups = group_structures['CASMO']['{}-group'.format(1)]
mgxs_lib = mgxs_lib.get_condensed_library(coarse_groups)

# Extract MGXS objects from library
nufiss = mgxs_lib.get_mgxs(mgxs_lib.domains[0], 'nu-fission')
total = mgxs_lib.get_mgxs(mgxs_lib.domains[0], 'total')
scatt = mgxs_lib.get_mgxs(mgxs_lib.domains[0], 'nu-scatter matrix')

# Get NumPy arrays of data
nufiss = nufiss.get_xs(nuclides='sum')
total = total.get_xs(nuclides='sum')
scatt = scatt.get_xs(nuclides='sum')
absorb = total - scatt

keff = nufiss / absorb
print('keff = {0:1.6f}'.format(float(keff[0])))
