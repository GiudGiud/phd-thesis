import os
import re

import matplotlib.pyplot as plt
import seaborn as sns

import pyne.ace

sns.set_style('ticks')

# Instantiate a PyNE ACE continuous-energy cross sections library
xs_dir = os.environ['OPENMC_CROSS_SECTIONS']
xs_dir = re.sub('cross_sections.xml$', '', xs_dir)
u238_lib = pyne.ace.Library(xs_dir + '293.6K/U_238_293.6K.ace')
u238_lib.read('92238.71c')

# Extract the U-238 data from the library
u238 = u238_lib.tables['92238.71c']

# Extract the continuous-energy cross section data
ce_capture = u238.reactions[102]

###########################  U-238 CAPTURE  ###############################
fig = plt.figure()
ax1 = fig.add_subplot(111)
plot1 = ax1.loglog(u238.energy, ce_capture.sigma,
                   color='mediumblue', linewidth=2)

# Begin plot customization
ax1.set_xlabel('Energy [MeV]', fontsize=12)
ax1.set_ylabel(r'$\sigma_{\gamma}$ [barns]', fontsize=12)

# Customize the plot
sns.set_style('ticks')
plt.title('U-238 Capture Cross Section', y=1.03, fontsize=16)
plt.xlim((1e-9, 1e1))
plt.savefig('u238-capture-xs.png', bbox_inches='tight')