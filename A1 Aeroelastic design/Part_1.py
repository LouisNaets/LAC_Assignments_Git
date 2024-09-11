import numpy as np
import math as m
import matplotlib.pyplot as plt

alpha_vary = np.arange(3, 7.1, 0.5) #[deg]
alpha = np.deg2rad(5) #[rad]

TSR_vary = np.arange(6, 10.1, 0.5) #[-]
TSR = 8 #[-]

B_vary = np.arange(2, 5.1, 1) #[-]
B = 3 #[-]

r_R = np.arange(0.1, 1.01, 0.01)

c_R_list = []
theta_R_list = []
for alpha_deg in alpha_vary:
    alpha_rad = np.deg2rad(alpha_deg)

    Cl_design = (2*m.pi)**2/360 * alpha_rad + 0.452

    c_R_list.append((16*m.pi/9) * (1 / (Cl_design * TSR**2 * B)) * 1/r_R)
    theta_R_list.append(2/3 * (1/TSR * 1/r_R) - alpha_rad) # [rad]

plt.style.use("seaborn-v0_8-whitegrid")
plt.rc('legend', fontsize=14)
plt.rc('axes', labelsize=16)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

fig, axs = plt.subplots(1, 2, figsize=(6.4, 7), dpi=100)

labels = ["%1.2f"%val for val in alpha_vary]

# Plot c/R vs r/R
for i, c_R in enumerate(c_R_list):
    axs[0].plot(r_R, c_R, label=f"α = {labels[i]}°")

axs[0].set_ylabel("c/R [-]")
axs[0].set_xlabel("r/R [-]")
axs[0].set_title("Chord distribution")
axs[0].legend()

# Plot θ/R vs r/R
for i, theta_R in enumerate(theta_R_list):
    axs[1].plot(r_R, np.rad2deg(theta_R), label=f"α = {labels[i]}°")

axs[1].set_ylabel("Twist [deg]")
axs[1].set_xlabel("r/R [-]")
axs[1].set_title("Twist distribution")
axs[1].legend()

fig.tight_layout()
plt.show()