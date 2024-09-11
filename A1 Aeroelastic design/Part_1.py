import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib.cm as cm

alpha_ref = 5 #[deg]
TSR_ref = 8 #[-]
B_ref = 3 #[-]

alpha_vary = np.arange(3, 7.1, 0.5) #[deg]
TSR_vary = np.arange(6, 10.1, 0.5) #[-]
B_vary = np.arange(2, 5.1, 1) #[-]

r_R = np.arange(0.1, 1.01, 0.01)

# Function to plot with color map and dotted lines
def plot_with_colormap(ax, x, y, label, cmap, idx, total_lines):
    colors = cm.viridis(np.linspace(0, 1, total_lines))  # Gradual colormap
    ax.plot(x, y, label=label, color=colors[idx], linestyle='-') 

fig, axs = plt.subplots(3, 2, figsize=(18, 12), dpi = 100)

for idx, alpha in enumerate(alpha_vary):
    r_R = np.arange(0.1, 1.01, 0.01)
    Cl_design = (2 * m.pi) ** 2 / 360 * alpha + 0.452
    c_R = (16 * m.pi / 9) * (1 / (Cl_design * TSR_ref ** 2 * B_ref)) * 1 / r_R
    theta_R = 2 / 3 * (1 / TSR_ref * 1 / r_R) - np.deg2rad(alpha)  # [rad]

    plot_with_colormap(axs[0, 0], r_R, c_R, label=f'α = {alpha}°', cmap=cm.viridis, idx=idx, total_lines=len(alpha_vary))
    plot_with_colormap(axs[0, 1], r_R, np.rad2deg(theta_R), label=f'α = {alpha}°', cmap=cm.viridis, idx=idx, total_lines=len(alpha_vary))

axs[0, 0].set_title("Chord distribution")
axs[0, 0].set_ylabel('c/R [-]')
axs[0, 0].grid(True, linestyle=':', linewidth=0.5)
axs[0, 0].legend()

axs[0, 1].set_title("Twist distribution")
axs[0, 1].set_ylabel('Twist [°]')
axs[0, 1].grid(True, linestyle=':', linewidth=0.5)
axs[0, 1].legend()

for idx, TSR in enumerate(TSR_vary):
    r_R = np.arange(0.1, 1.01, 0.01)
    Cl_design = (2 * m.pi) ** 2 / 360 * alpha_ref + 0.452
    c_R = (16 * m.pi / 9) * (1 / (Cl_design * TSR ** 2 * B_ref)) * 1 / r_R
    theta_R = 2 / 3 * (1 / TSR * 1 / r_R) - np.deg2rad(alpha_ref)  # [rad]

    plot_with_colormap(axs[1, 0], r_R, c_R, label=f'λ = {TSR}', cmap=cm.viridis, idx=idx, total_lines=len(TSR_vary))
    plot_with_colormap(axs[1, 1], r_R, np.rad2deg(theta_R), label=f'λ = {TSR}', cmap=cm.viridis, idx=idx, total_lines=len(TSR_vary))

axs[1, 0].set_ylabel('c/R [-]')
axs[1, 0].grid(True, linestyle=':', linewidth=0.5)
axs[1, 0].legend()

axs[1, 1].set_ylabel('Twist [°]')
axs[1, 1].grid(True, linestyle=':', linewidth=0.5)
axs[1, 1].legend()

for idx, B in enumerate(B_vary):
    r_R = np.arange(0.1, 1.01, 0.01)
    Cl_design = (2 * m.pi) ** 2 / 360 * alpha_ref + 0.452
    c_R = (16 * m.pi / 9) * (1 / (Cl_design * TSR_ref ** 2 * B)) * 1 / r_R
    theta_R = 2 / 3 * (1 / TSR_ref * 1 / r_R) - np.deg2rad(alpha_ref)  # [rad]

    plot_with_colormap(axs[2, 0], r_R, c_R, label=f'B = {B}', cmap=cm.viridis, idx=idx, total_lines=len(B_vary))
    plot_with_colormap(axs[2, 1], r_R, np.rad2deg(theta_R), label=f'B = {B}', cmap=cm.viridis, idx=idx, total_lines=len(B_vary))

axs[2, 0].set_xlabel('r/R [-]')
axs[2, 0].set_ylabel('c/R [-]')
axs[2, 0].grid(True, linestyle=':', linewidth=0.5)
axs[2, 0].legend()

axs[2, 1].set_xlabel('r/R [-]')
axs[2, 1].set_ylabel('Twist [°]')
axs[2, 1].grid(True, linestyle=':', linewidth=0.5)
axs[2, 1].legend()

fig.tight_layout()
plt.show()

phi = np.rad2deg(2/3*(1/TSR_ref*1/r_R))

fig2 = plt.figure(1, figsize=(6, 4), dpi = 100)
plt2 = fig2.add_subplot()
plt2.plot(r_R, phi)
plt2.grid()

#plt2.set_title('Inflow angle')
plt2.set_xlabel('r/R [-]')
plt2.set_ylabel('Inflow angle [°]')

fig2.tight_layout()
plt.show()