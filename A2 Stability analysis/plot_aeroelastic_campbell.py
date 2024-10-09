"""Plot an aeroelastic Campbell diagram versus wind speed.
"""
import matplotlib.pyplot as plt
import numpy as np

from lacbox.io import load_cmb, load_amp
from lacbox.vis import plot_amp
from lacbox.io import load_st


TURBINE_NAME = 'Redesigned IIIB climate turbine'
CMB_PATH = './A2/group7_aeroelastic.cmb'
NMODES = 8  # number of modes to plot
MODE_NAMES = ['Tower side-side', 'Tower fore-aft', '1st flap BW', '1st flap FW', '1st flap SYM',
              '1st edge BW', '1st edge FW', '2nd flap BW', '2nd flap FW', '2nd flap SYM',
              '1st edge SYM']
OPT_PATH = None  # path to opt file, needed for P-harmonics

# load campbell diagram
wsp, dfreqs, zetas = load_cmb(CMB_PATH, cmb_type='aeroelastic')

# Define custom color and marker groups for specific mode ranges
color_marker_map = {
    (0, 1): ('tab:blue', ['o', 'x']), 
    (2, 3, 4): ('tab:orange', ['o', 'x', '^']),  
    (5, 6, 10): ('tab:green', ['o', 'x', '^']),   
    (7, 8, 9): ('tab:red', ['o', 'x', '^'])             
}

# initialize plot
fig, axs = plt.subplots(1, 2, figsize=(9.5, 4), dpi=500)

# loop through modes
NMODES = len(MODE_NAMES)
for i in range(NMODES):
    # Set default marker and color
    m = '.'
    c = 'tab:gray'
    # Find the corresponding color and marker for this mode
    for mode_range, (color, markers) in color_marker_map.items():
        if i in mode_range:
            m = markers[mode_range.index(i)]  # Get the specific marker
            c = color  # Use the defined color
            break

    # left plot: damped nat freqs in ground-fixed frame
    axs[0].plot(wsp, dfreqs[:, i], marker=m, c=c, mfc='none', ms=5, linewidth=0.75, markeredgewidth=0.65)

    # right plot: percent criticl damping
    axs[1].plot(wsp, zetas[:, i], marker=m, label=MODE_NAMES[i], c=c, mfc='none', ms=5, linewidth=0.75, markeredgewidth=0.65)

# load opt file, add P-harmonics?
opt_path = './A2/group7_3B_design_flex.opt'
flex_opt_data = np.loadtxt(opt_path, skiprows=1)

class f_opt:
    V_o = flex_opt_data[:, 0]
    omega = flex_opt_data[:, 2]

#Plotting P harmonics 1, 3, and 6
P_harmonics = [f_opt.omega/60*1, f_opt.omega/60*3, f_opt.omega/60*6]
for i in range(len(P_harmonics)):
    axs[0].plot(f_opt.V_o, P_harmonics[i], linewidth=0.85, c='tab:gray')
    if i == 1:
        axs[0].plot(f_opt.V_o, P_harmonics[i], linewidth=0.85, c='tab:gray', label='1P, 3P, 6P')
        axs[1].plot(5, 0, linewidth=0.85, c='tab:gray', label='1P, 3P, 6P') #for the legend

# prettify
axs[0].set(xlabel='Wind speed [m/s]', ylabel='Damped nat. frequencies [Hz]')
axs[0].grid()
axs[1].set(xlabel='Wind speed [m/s]', ylabel='Modal damping [% critical]')
axs[1].grid()
axs[1].legend(bbox_to_anchor=(1.02, 0.5), loc='center left')

# add figure title and scale nicely
fig.suptitle(f'Aeroelastic Campbell diagram for {TURBINE_NAME}')
fig.tight_layout()

fig.savefig('A2/Figures/Campbell_aeroelastic.svg', format='svg')
fig.savefig('A2/Figures/Campbell_aeroelastic.png', format='png')

nmodes = dfreqs.shape[1]  # get number of modes
mode_names = [f'Mode {i}' for i in range(1, nmodes+1)]  # list of mode shape names

# Path to the .amp file
amp_path = './A2/group7_aeroelastic_amp.amp'

# Load the modal amplitudes
amp_df = load_amp(amp_path)
#print(amp_df.shape)
#print(amp_df)

#print(amp_df.index)
#amp_df.head()

wsp = 11.1  # rated wind speed

fig, ax = plot_amp(amp_df, mode_names, wsp, title=f'Redesigned IIIB climate turbine aeroelastic modal amplitudes at ~{wsp:.0f} m/s')

fig.savefig('A2/Figures/Campbell_aeroelastic_modal_amplitudes.svg', format='svg')
fig.savefig('A2/Figures/Campbell_aeroelastic_modal_amplitudes.png', format='png')