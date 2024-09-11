import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Parameters
aoa_design_ref = np.deg2rad(5)  # [rad]
B_ref = 3  # [-]
TSR_ref = 8  # [-]

aoa_design_range = np.deg2rad(np.arange(3, 8, 1))
B_range = np.arange(2, 6, 1)
TSR_range = np.arange(6, 11, 1)

# Function to plot with color map and dotted lines
def plot_with_colormap(ax, x, y, label, cmap, idx, total_lines):
    colors = cm.viridis(np.linspace(0, 1, total_lines))  # Gradual colormap
    ax.plot(x, y, label=label, color=colors[idx], linestyle='-') 

fig, axes = plt.subplots(3, 2, figsize=(18, 12))

# Plot for TSR range (Chord Distribution and Twist Distribution)
for idx, TSR in enumerate(TSR_range):
    r_R = np.arange(0.1, 1.01, 0.01)
    Cl_design = (2 * m.pi) ** 2 / 360 * aoa_design_ref + 0.452
    c_R = (16 * m.pi / 9) * (1 / (Cl_design * TSR ** 2 * B_ref)) * 1 / r_R
    theta_R = 2 / 3 * (1 / TSR * 1 / r_R) - aoa_design_ref  # [rad]

    plot_with_colormap(axes[0, 0], r_R, c_R, label=f'TSR = {TSR}', cmap=cm.viridis, idx=idx, total_lines=len(TSR_range))
    plot_with_colormap(axes[0, 1], r_R, np.rad2deg(theta_R), label=f'TSR = {TSR}', cmap=cm.viridis, idx=idx, total_lines=len(TSR_range))

# axes[0, 0].set_xlabel('r/R [-]')
axes[0, 0].set_ylabel('c/R [-]')
axes[0, 0].set_title('Chord Distribution')
axes[0, 0].grid(True, linestyle=':', linewidth=0.5)
axes[0, 0].legend()

# axes[0, 1].set_xlabel('r/R [-]')
axes[0, 1].set_ylabel('Twist [°]')
axes[0, 1].set_title('Twist Distribution')
axes[0, 1].grid(True, linestyle=':', linewidth=0.5)
axes[0, 1].legend()

# Plot for AoA range (Chord Distribution and Twist Distribution)
for idx, aoa in enumerate(aoa_design_range):
    r_R = np.arange(0.1, 1.01, 0.01)
    Cl_design = (2 * m.pi) ** 2 / 360 * aoa + 0.452
    c_R = (16 * m.pi / 9) * (1 / (Cl_design * TSR_ref ** 2 * B_ref)) * 1 / r_R
    theta_R = 2 / 3 * (1 / TSR_ref * 1 / r_R) - aoa  # [rad]

    plot_with_colormap(axes[1, 0], r_R, c_R, label=f'aoa =  {int(np.rad2deg(aoa))}°', cmap=cm.plasma, idx=idx, total_lines=len(aoa_design_range))
    plot_with_colormap(axes[1, 1], r_R, np.rad2deg(theta_R), label=f'aoa = {int(np.rad2deg(aoa))}°', cmap=cm.plasma, idx=idx, total_lines=len(aoa_design_range))

# axes[1, 0].set_xlabel('r/R [-]')
axes[1, 0].set_ylabel('c/R [-]')
axes[1, 0].grid(True, linestyle=':', linewidth=0.5)
axes[1, 0].legend()

# axes[1, 1].set_xlabel('r/R [-]')
axes[1, 1].set_ylabel('Twist [°]')
axes[1, 1].grid(True, linestyle=':', linewidth=0.5)
axes[1, 1].legend()

# Plot for B range (Chord Distribution and Twist Distribution)
for idx, B in enumerate(B_range):
    r_R = np.arange(0.1, 1.01, 0.01)
    Cl_design = (2 * m.pi) ** 2 / 360 * aoa_design_ref + 0.452
    c_R = (16 * m.pi / 9) * (1 / (Cl_design * TSR_ref ** 2 * B)) * 1 / r_R
    theta_R = 2 / 3 * (1 / TSR_ref * 1 / r_R) - aoa_design_ref  # [rad]

    plot_with_colormap(axes[2, 0], r_R, c_R, label=f'B = {B}', cmap=cm.cividis, idx=idx, total_lines=len(B_range))
    plot_with_colormap(axes[2, 1], r_R, np.rad2deg(theta_R), label=f'B = {B}', cmap=cm.cividis, idx=idx, total_lines=len(B_range))

axes[2, 0].set_xlabel('r/R [-]')
axes[2, 0].set_ylabel('c/R [-]')
axes[2, 0].grid(True, linestyle=':', linewidth=0.5)
axes[2, 0].legend()

axes[2, 1].set_xlabel('r/R [-]')
axes[2, 1].set_ylabel('Twist [°]')
axes[2, 1].grid(True, linestyle=':', linewidth=0.5)
axes[2, 1].legend()

plt.tight_layout()
plt.show()
