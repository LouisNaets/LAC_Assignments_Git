"""Plot a structural Campbell diagram versus rotor speed.

Note:
    - The special variable to indicate which index is the drivetrain mode. Give
      the mode number according to _HAWCStab2_.
    - Structural Campbell diagrams often incorrectly label modes after crossings.
      I added logic to "swap" the modes back after a certain RPM. You may need to
      update that logic for your turbine. Verify with HAWCStab2 animations.
"""
import matplotlib.pyplot as plt
import numpy as np

from lacbox.io import load_cmb, load_amp
from lacbox.vis import plot_amp
from plot_amp_custom import plot_amp_rpm


TURBINE_NAME = 'Redesigned IIIB climate turbine'
CMB_PATH = './A2/group7_structural.cmb'
MODE_NAMES = ['Tower side-side', 'Tower fore-aft', '1st flap BW', '1st flap FW', '1st flap SYM',
              '1st edge BW', '1st edge FW', '2nd flap BW', '2nd flap FW', '2nd flap SYM',
              '1st edge SYM']
DT_MODENUM = 9  # what mode number in HAWCStab2 is the drivetrain mode?
# define turbine natural frequency? blade-only natural frequencies?

# load campbell diagram
omega, dfreqs, _ = load_cmb(CMB_PATH, cmb_type='structural')
omega_rpm = omega * 30 / np.pi  # from rad/s to RPM

# swap modes 5 and 6 after 3 RPM. hawcstab2 mixed them up at the crossing.
# note that numpy index is off by 2 because of python + skipping rigid-body mode.
SWAP_RPM = 3
dfreq_mode5 = dfreqs[omega_rpm > SWAP_RPM, 3]
dfreqs[omega_rpm > SWAP_RPM, 3] = dfreqs[omega_rpm > SWAP_RPM, 4]
dfreqs[omega_rpm > SWAP_RPM, 4] = dfreq_mode5

# Define custom color and marker groups for specific mode ranges
color_marker_map = {
    (0, 1): ('tab:blue', ['o', 'x']), 
    (2, 3, 4): ('tab:orange', ['o', 'x', '^']),  
    (5, 6, 10): ('tab:green', ['o', 'x', '^']),   
    (7, 8, 9): ('tab:red', ['o', 'x', '^'])             
}

# Initialize plot
fig, ax = plt.subplots(figsize=(8, 4), dpi=500)

# Loop through modes
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

    # If we're at the end of the loop, use the drivetrain mode number index
    arr_idx = i if i < NMODES - 1 else DT_MODENUM - 2  # Off by 2 from HAWCStab2 due to skipping rigid-body

    # Plot damped natural frequencies in ground-fixed frame
    ax.plot(omega_rpm, dfreqs[:, arr_idx], marker=m, c=c, label=MODE_NAMES[i], mfc='none', ms=5, linewidth=0.75, markeredgewidth=0.65)

#New color marker map
color_marker_map = {
    (0, 1): ('tab:blue', ['o', 'x']), 
    (2, 3, 4): ('tab:orange', ['o', 'x', '^']),  
    (5, 6, 7): ('tab:green', ['o', 'x', '^'])           
}

#Plotting the simple versions
natural_freqs_1 = [[(0.25+0.25)/2], [(0.52 + 0.57 + 0.61)/3], [(0.88 + 0.89)/2]]
natural_freqs = []
for i in range(0,len(natural_freqs_1)):
  natural_freqs.append(np.repeat(natural_freqs_1[i], len(omega_rpm)))

omega_hz=omega_rpm/60

simple_freqs = [natural_freqs[0], natural_freqs[0], natural_freqs[1]-omega_hz, natural_freqs[1]+omega_hz, natural_freqs[1], natural_freqs[2]-omega_hz, natural_freqs[2]+omega_hz, natural_freqs[2]]

for i in range(len(simple_freqs)):
    # Set default marker and color
    m = '.'
    c = 'tab:gray'
    # Find the corresponding color and marker for this mode
    for mode_range, (color, markers) in color_marker_map.items():
        if i in mode_range:
            m = markers[mode_range.index(i)]  # Get the specific marker
            c = color  # Use the defined color
            break
    ax.plot(omega_rpm, simple_freqs[i], c=c, linewidth=0.75, linestyle='--')#, mfc='none', ms=5, marker=m, markeredgewidth=0.65)
ax.plot(0, 0, c='black', linewidth=0.75, linestyle='--',label = 'Theory')

#Plotting P harmonics 1, 3, and 6
P_harmonics = [omega_rpm/60*1, omega_rpm/60*3, omega_rpm/60*6]
for i in range(len(P_harmonics)):
    ax.plot(omega_rpm, P_harmonics[i], linewidth=0.85, c='tab:gray')
    if i == 1:
        ax.plot(omega_rpm, P_harmonics[i], linewidth=0.75, c='tab:gray', label='1P, 3P, 6P')

# prettify
ax.set(xlabel='Rotor speed [RPM]', ylabel='Damped nat. frequencies [Hz]')
ax.grid()
ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left')

# add figure title and scale nicely
fig.suptitle(f'Structural Campbell diagram for {TURBINE_NAME}')
fig.tight_layout()

fig.savefig('A2/Figures/Campbell_structural.svg', format='svg')
fig.savefig('A2/Figures/Campbell_structural.png', format='png')

nmodes = dfreqs.shape[1]  # get number of modes
mode_names = [f'Mode {i}' for i in range(1, nmodes+1)]  # list of mode shape names

# Path to the .amp file
amp_path = './A2/group7_structural_amp.amp'

# Load the modal amplitudes
amp_df = load_amp(amp_path)

#print(amp_df.index)
#amp_df.head()

rpm = 0.279253E+00  # rotational speed

fig, ax = plot_amp_rpm(amp_df, mode_names, rpm, title=f'Redesigned IIIB climate turbine structural modal amplitudes at ~{rpm*60/(2*np.pi):.0f} RPM')

fig.savefig('A2/Figures/Campbell_structural_modal_amplitudes.svg', format='svg')
fig.savefig('A2/Figures/Campbell_structural_modal_amplitudes.png', format='png')