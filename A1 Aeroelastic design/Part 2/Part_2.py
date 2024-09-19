# %% Import modules
import matplotlib.pyplot as plt
import numpy as np
from aero_design_functions import get_design_functions, single_point_design
from sympy import symbols, Eq, solve

"""
DTU 10MW Class IA -> IIIB
"""
# Global constants
rho = 1.225 # Air Density [kg/m^3]
P_rated = 10e6 # Rated Power [W]
tsr = 9.0 # Tip-Speed-Ratio [-]
B = 3  # Number of blades [#]
cl_scale = 1 # Change this value to scale the cl-values

# IEC 61400-1 Classes
class IA:
    I = 0.18 # [-]
    V_rated = 11.4 # m/s 
    R = 89.15 # m
    r_hub = 2.8  # Hub radius [m]
    chord_max = 6.20  # Maximum chord size [m]
    chord_root = 5.38  # Chord size at the root [m]

class IIIB:
    I = 0.16 # [-]


#%% Step 2: Rotor radius and rated wind speed for IIIB class
# System of equations solver for IIIB class R and V_rated
IIIB_R, IIIB_V_rated = symbols('IIIB_R IIIB_V_rated')
eq1 = ((IA.V_rated * (1 + 2 * IA.I)) / (IIIB_V_rated * (1 + 2*IIIB.I)))**(2/3) * IA.R - IIIB_R
eq2 = (IA.R / IIIB_R)**(2/3) * IA.V_rated - IIIB_V_rated
solution = solve((eq1, eq2), (IIIB_R, IIIB_V_rated))[0]

# Assign the solutions to the IIIB object
IIIB.R, IIIB.V_rated = solution
IIIB.R = float(IIIB.R)
IIIB.V_rated = float(IIIB.V_rated)

print(f"Rotor radius for IIIB class: {IIIB.R:.2f} m")
print(f"Rated wind speed for IIIB class: {IIIB.V_rated:.2f} m/s")

#%% Step 4: Upscaling the blade thickness and centerline chord length

# Inputs scaled to IIIB class 10 MW turbine
SF = IIIB.R / IA.R # Scaling factor Assuming linear scaling
IIIB.r_hub = SF * IA.r_hub  # Hub radius [m]
IIIB.r = np.linspace(IIIB.r_hub, IIIB.R - 0.1, 40)  # Rotor span [m]
IIIB.chord_max = SF * IA.chord_max  # Maximum chord size [m]
IIIB.chord_root = SF * IA.chord_root  # Chord size at the root [m]

# Determine chord and thickness for IA class from the _ae file
file_path = r'hawc_files\dtu_10mw\data\DTU_10MW_RWT_ae.dat'
data = np.loadtxt(file_path, skiprows=2, usecols=(0, 1, 2))
IA.r = data[:, 0]  # Blade span [m]
IA.chord = data[:, 1]  # Chord length [m]
IA.tc = data[:, 2]/100  # Thickness-to-chord ratio [-]
IA.t = IA.tc * IA.chord  # Absolute Thickness [m]

# Upscale to IIIB class
IIIB.r = SF * IA.r
IIIB.chord = SF * IA.chord
IIIB.tc = SF * IA.tc
IIIB.t = SF*IA.t # Absolute Thickness [m]

print(f"Chord size at the root for IIIB class: {IIIB.chord_root:.2f} m")
print(f"Maximum chord size for IIIB class: {IIIB.chord_max:.2f} m")

# Plot absolute thickness distributions
fig2, ax2 = plt.subplots(1, 1, num=2)
ax2.plot(IA.r, IA.t, label='DTU 10MW Class IA')
ax2.plot(IIIB.r, IIIB.t, label='DTU 10MW Class IIIB')
ax2.set_xlabel('Blade Span [m]')
ax2.set_ylabel('Absolute Thickness [m]')
ax2.legend()
fig2.tight_layout()

#%% Step 6: Choice of design lift and airfoils for 3 different airfoils

pc_path = r'hawc_files\dtu_10mw\data\DTU_10MW_RWT_pc.dat'
tc_data = np.loadtxt(pc_path, skiprows=2, usecols=[0,1,2])
tc_241_data = tc_data[52-2:69-2, :]
tc_301_data = tc_data[158-2:175-2, :]
tc_360_data = tc_data[264-2:281-2, :]
tc_480_data = tc_data[370-2:379-2, :]

# List of tc datasets for convenience
tc_datasets = [tc_241_data, tc_301_data, tc_360_data, tc_480_data]
titles = ['tc_241_data', 'tc_301_data', 'tc_360_data', 'tc_480_data']

# Plot for each tc dataset
# for i, tc_data in enumerate(tc_datasets):
#     alpha = tc_data[:, 0]  # Alpha (angle of attack)
#     c_l = tc_data[:, 1]    # Lift coefficient (C_l)
#     c_d = tc_data[:, 2]    # Drag coefficient (C_d)
#     c_l_c_d = c_l / c_d    # Lift-to-drag ratio (C_l/C_d)

#     # Create a figure for each tc_data set
#     fig, axs = plt.subplots(1, 3, figsize=(18, 6))
#     fig.suptitle(f"Plots for {titles[i]}")
    
#     # Plot 1: C_l vs C_d
#     axs[0].plot(c_d, c_l, 'o-', color='blue')
#     axs[0].set_title(f'C_l vs C_d')
#     axs[0].set_xlabel('C_d')
#     axs[0].set_ylabel('C_l')
#     axs[0].grid()

#     # Plot 2: C_l vs Alpha
#     axs[1].plot(alpha, c_l, 'o-', color='green')
#     axs[1].set_title(f'C_l vs Alpha')
#     axs[1].set_xlabel('Alpha [deg]')
#     axs[1].set_ylabel('C_l')
#     axs[1].grid()

#     # Plot 3: C_l vs C_l/C_d
#     axs[2].plot(c_l_c_d, c_l, 'o-', color='red')
#     axs[2].set_title(f'C_l vs C_l/C_d')
#     axs[2].set_xlabel('C_l/C_d')
#     axs[2].set_ylabel('C_l')
#     axs[2].grid()

#     # Display the plots
#     plt.tight_layout()

# Chosen design points for each airfoil
Cl_points = [1.3, 1.23, 1.37, 0.54]
Cd_points = [0.013, 0.014, 0.021, 0.032]
alpha_points = [8.1, 7.5, 5.8, 1.8]
t_c_points = [24.1, 30.1, 36.0, 48.0]
design_points = dict(zip(titles, zip(Cl_points, alpha_points)))

#%%

# DTU 10MW IA
fig1, axs1 = plt.subplots(3, 1, num=1)
fig2, axs2 = plt.subplots(3, 1, num=2, clear=True)
fig3, axs3 = plt.subplots(2, 2, num=3, clear=True)
fig4, axs4 = plt.subplots(3, 1, num=4, clear=True, figsize=(6.5, 5.5))

i = np.arange(1, 4)
colors = ['tab:blue', 'tab:orange', 'tab:green']
for i_design in i:
    IA.r = np.linspace(IA.r_hub, IA.R - 0.1, 40)  # Rotor span [m]
    IA.cl_des, IA.cd_des, IA.aoa_des, IA.tc_vals, IA.cl_vals, IA.cd_vals, IA.aoa_vals = get_design_functions(i_design)

    IA.chord, IA.tc, IA.twist, IA.cl, IA.cd, IA.aoa, IA.a, IA.CLT, IA.CLP, IA.CT, IA.CP = single_point_design(
        IA.r, IA.t, tsr, IA.R, IA.cl_des, IA.cd_des, IA.aoa_des, IA.chord_root, IA.chord_max, B
    )

    # %% Plotting design functions
    tc_plot = np.linspace(0, 100, 100)

    axs1[0].plot(tc_plot, IA.cl_des(tc_plot), color = colors[i_design-1], label=f'Des_i {i_design}')
    axs1[0].plot(IA.tc_vals, IA.cl_vals, "o", color = colors[i_design-1])
    axs1[0].plot(t_c_points, Cl_points, 'ok')
    axs1[0].set_ylabel("$C_l$ [-]")
    axs1[0].set_xlim(0, 100)
    axs1[0].legend()

    axs1[1].plot(tc_plot, IA.cd_des(tc_plot), color = colors[i_design-1], label=f'Des_i {i_design}')
    axs1[1].plot(IA.tc_vals, IA.cd_vals, "o", color = colors[i_design-1])
    axs1[1].plot(t_c_points, Cd_points, 'ok')
    axs1[1].set_ylabel("$C_d$ [-]")
    axs1[1].set_xlim(0, 100)
    axs1[1].legend()

    axs1[2].plot(tc_plot, IA.aoa_des(tc_plot), color = colors[i_design-1], label=f'Des_i {i_design}')
    axs1[2].plot(IA.tc_vals, IA.aoa_vals, "o", color = colors[i_design-1])
    axs1[2].plot(t_c_points, alpha_points, 'ok')
    axs1[2].set_ylabel(r"$\alpha$ [-]")
    axs1[2].set_xlabel(r"$t/c$ [-]")
    axs1[2].set_xlim(0, 100)
    axs1[2].legend()

    fig1.tight_layout()

    # %% Plot the chord, twist and relative-thickness
    # Chord
    axs2[0].plot(IA.r, IA.chord, label=f'Des_i {i_design}')
    axs2[0].set_ylabel("Chord [m]")
    axs2[0].set_xlim(0, max(IA.R, IA.R))
    axs2[0].legend()

    # Twist
    axs2[1].plot(IA.r, IA.twist, label=f'Des_i {i_design}')
    axs2[1].set_ylabel("Twist [deg]")
    axs2[1].set_xlim(0, max(IA.R, IA.R))
    axs2[1].legend()

    # t/c
    axs2[2].plot(IA.r, IA.tc, label=f'Des_i {i_design}')
    axs2[2].set_ylabel("Rel. thickness [%]")
    axs2[2].set_xlabel("Rotor span [m]")
    axs2[2].set_xlim(0, max(IA.R, IA.R))
    axs2[2].legend()

    # %% Plot r vs. t/c, aoa, cl, cd
    # t/c
    axs3[0, 0].plot(IA.r, IA.tc, label=f'Des_i {i_design}')
    axs3[0, 0].set_ylabel("t/c [%]")
    axs3[0, 0].set_xlim(0, max(IA.R, IA.R))
    axs3[0, 0].legend()

    # aoa
    axs3[0, 1].plot(IA.r, IA.aoa, label=f'Des_i {i_design}')
    axs3[0, 1].set_ylabel(r"$\alpha$ [deg]")
    axs3[0, 1].set_xlim(0, max(IA.R, IA.R))
    axs3[0, 1].yaxis.tick_right()
    axs3[0, 1].yaxis.set_label_position("right")
    axs3[0, 1].legend()

    # cl
    axs3[1, 0].plot(IA.r, IA.cl, label=f'Des_i {i_design}')
    axs3[1, 0].set_ylabel("$C_l$ [-]")
    axs3[1, 0].set_xlabel("Span [m]")
    axs3[1, 0].set_xlim(0, max(IA.R, IA.R))
    axs3[1, 0].legend()

    # cd
    axs3[1, 1].plot(IA.r, IA.cd, label=f'Des_i {i_design}')
    axs3[1, 1].set_ylabel("$C_d$ [-]")
    axs3[1, 1].set_xlabel("Span [m]")
    axs3[1, 1].set_xlim(0, max(IA.R, IA.R))
    axs3[1, 1].yaxis.tick_right()
    axs3[1, 1].yaxis.set_label_position("right")
    axs3[1, 1].legend()

    # %% Plot r vs. CLT, CLP, a
    # Local-Thrust-Coefficient
    axs4[0].plot(IA.r, IA.CLT, label=f'Des_i {i_design}')
    axs4[0].axhline(y=8 / 9, ls="--", color="k", lw=1)
    axs4[0].set_ylabel("Local thrust ($C_{LT}$) [-]")
    axs4[0].set_ylim(0, 1.0)
    axs4[0].set_xlim(0, max(IA.R, IA.R))
    axs4[0].legend()

    # Local-Power-Coefficient
    axs4[1].plot(IA.r, IA.CLP, label=f'Des_i {i_design}')
    axs4[1].axhline(y=16 / 27, ls="--", color="k", lw=1)
    axs4[1].set_ylabel("Local Power ($C_{LP}$) [-]")
    axs4[1].set_xlim(0, max(IA.R, IA.R))
    axs4[1].set_ylim(-0.4, 0.6)
    axs4[1].legend()

    # Axial Induction
    axs4[2].plot(IA.r, IA.a, label=f'Des_i {i_design}')
    axs4[2].axhline(y=1 / 3, ls="--", color="k", lw=1)
    axs4[2].set_ylabel("Axial induction ($a$) [-]")
    axs4[2].set_xlabel("Rotor span [m]")
    axs4[2].set_xlim(0, max(IA.R, IA.R))
    axs4[2].legend()

    fig4.suptitle(f"IA: $C_T$={IA.CT:1.3f}, $C_P$={IA.CP:1.3f}")
    
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()

plt.show()
