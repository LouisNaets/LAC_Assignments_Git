"""Plot spanwise multiple TSRs with a nice colormap.
"""
from pathlib import Path

from lacbox.io import load_ind
import matplotlib as mpl
import matplotlib.pyplot as plt


IND_DIR = (r'C:\Users\rink\git\G-teaching\lac-course-private\dtu-10mw'
           + r'\res_hawc2s_multitsr')  # directory with induction files
FILENAME_BASE = 'dtu_10mw_hawc2s_multitsr_u'  # beginning of the .ind files we want
SKIP_FILES = 3  # I calculated for a lot of TSRs, so I want to skip some of them for a clean plot
XLIM = [0, 87]  # hard-code my x-limits
SAVE_FIG = True  # save the figure as a png?

# get sorted list of induction files, skipping some if needed
ind_genr = Path(IND_DIR).glob(FILENAME_BASE + '*')  # returns a generator
ind_paths = sorted(list(ind_genr))[::SKIP_FILES]  # convert to list and sort (just in case)
n_lines = len(ind_paths)

# initialize the figure
fig, axs = plt.subplots(3, 2, figsize=(6.4, 6))

# choose our colormap
cmap = plt.get_cmap('viridis')

# loop over induction files
for i, ind_path in enumerate(ind_paths):

    # choose the color for this file
    line_color = cmap(i/(n_lines - 1))

    # load the data
    ind_data = load_ind(ind_path)

    # a
    axs[0, 0].plot(ind_data["s_m"], ind_data["a"], c=line_color)
    axs[0, 0].set(ylabel="ax. ind. ($a$) [-]")

    # ap
    axs[0, 1].plot(ind_data["s_m"], ind_data["ap"], c=line_color)
    axs[0, 1].set(ylabel="tan. ind. ($a_p$) [-]")

    # Cl
    axs[1, 0].plot(ind_data["s_m"], ind_data["Cl"], c=line_color)
    axs[1, 0].set(ylabel="$C_l$ [-]")

    # Cd
    axs[1, 1].plot(ind_data["s_m"], ind_data["Cd"], c=line_color)
    axs[1, 1].set(ylabel="$C_d$ [-]")

    # CP
    axs[2, 0].plot(ind_data["s_m"], ind_data["CP"], c=line_color)
    axs[2, 0].set(xlabel="Blade-span ($s$) [m]",
                  ylabel="local $C_P$ [-]")

    # CP
    axs[2, 1].plot(ind_data["s_m"], ind_data["CT"], c=line_color)
    axs[2, 1].set(xlabel="Blade-span ($s$) [m]",
                  ylabel="local-$C_T$ [-]")

# set "global" values on axes
for ax in fig.get_axes():
    ax.grid()
    ax.set(xlim=XLIM)

# plt.colorbar(location='top', orientation='horizontal')
cax = fig.add_axes([0.15, 0.96, 0.7, 0.02])
norm = mpl.colors.Normalize(vmin=6, vmax=9)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=cax, orientation='horizontal', label='TSR')

fig.tight_layout(rect=[0, 0, 1, 0.9])

if SAVE_FIG:
    fig.savefig('multitsr.png', dpi=150)

plt.show()
