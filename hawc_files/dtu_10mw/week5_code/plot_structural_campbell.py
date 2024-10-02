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

from lacbox.io import load_cmb


TURBINE_NAME = 'DTU 10 MW'
CMB_PATH = './dtu_10mw_structural.cmb'
MODE_NAMES = ['Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', 'Mode 6',
              'Mode 7', 'Mode 8', 'Mode 9']
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

# initialize plot
fig, ax = plt.subplots(figsize=(6, 4))

# plot theoretical lines?

# loop through modes
NMODES = len(MODE_NAMES)
for i in range(NMODES):
    # add a custom color or marker for each mode?
    m, c = '.', f'C{i}'

    # if we're at the end of the loop, take the index of the drivetrain mode number instead
    arr_idx = i if i < NMODES - 1 else DT_MODENUM - 2  # off by 2 from hawcstab2 because python + skip rigid-body

    # damped nat freqs in ground-fixed frame
    ax.plot(omega_rpm, dfreqs[:, arr_idx], marker=m, c=c, label=MODE_NAMES[i])

# prettify
ax.set(xlabel='Rotor speed [RPM]', ylabel='Damped nat. frequencies [Hz]')
ax.grid()
ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left')

# add figure title and scale nicely
fig.suptitle(f'Structural Campbell diagram for {TURBINE_NAME}')
fig.tight_layout()

plt.show()
