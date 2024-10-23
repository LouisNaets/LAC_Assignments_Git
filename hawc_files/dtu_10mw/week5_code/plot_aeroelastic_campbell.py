"""Plot an aeroelastic Campbell diagram versus wind speed.
"""
import matplotlib.pyplot as plt
import numpy as np

from lacbox.io import load_cmb


TURBINE_NAME = 'DTU 10 MW'
CMB_PATH = 'dtu_10mw_aeroelastic.cmb'
NMODES = 8  # number of modes to plot
MODE_NAMES = ['Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', 'Mode 6',
              'Mode 7', 'Mode 8', 'Mode 9', 'Mode 10', 'Mode 11',
              'Mode 12']
OPT_PATH = None  # path to opt file, needed for P-harmonics

# load campbell diagram
wsp, dfreqs, zetas = load_cmb(CMB_PATH, cmb_type='aeroelastic')

# initialize plot
fig, axs = plt.subplots(1, 2, figsize=(9.5, 4))

# loop through modes
NMODES = len(MODE_NAMES)
for i in range(NMODES):
    # add a custom color or marker for each mode?
    m, c = '.', f'C{i}'

    # left plot: damped nat freqs in ground-fixed frame
    axs[0].plot(wsp, dfreqs[:, i], marker=m, c=c)

    # right plot: percent criticl damping
    axs[1].plot(wsp, zetas[:, i], marker=m, c=c, label=MODE_NAMES[i])

# load opt file, add P-harmonics?

# prettify
axs[0].set(xlabel='Wind speed [m/s]', ylabel='Damped nat. frequencies [Hz]')
axs[0].grid()
axs[1].set(xlabel='Wind speed [m/s]', ylabel='Modal damping [% critical]')
axs[1].grid()
axs[1].legend(bbox_to_anchor=(1.02, 0.5), loc='center left')

# add figure title and scale nicely
fig.suptitle(f'Aeroelastic Campbell diagram for {TURBINE_NAME}')
fig.tight_layout()
fig.savefig(f'./Figures/Aeroelastic_campbell.svg', format='svg')
plt.show()
