# %% Import modules
import matplotlib.pyplot as plt
import numpy as np
from aero_design_functions import get_design_functions, single_point_design
from sympy import symbols, Eq, solve
from matplotlib import cm

"""
Assignment 1: Aeroelastic Design part 2
DTU 10MW Class IA -> IIIB
"""

# Global constants
rho = 1.225 # Air Density [kg/m^3]
P_rated = 10e6 # Rated Power [W]
tsr = 9.0 # Tip-Speed-Ratio [-]
B = 3  # Number of blades [#]

# IEC 61400-1 Classes
class IA:
    I = 0.18 # Turbulence Intensity[-]
    V_rated = 11.4 # m/s 
    r_hub = 2.8  # Hub radius [m]
    R = 89.15 - r_hub # Rotor diameter (excl hub) [m]
    chord_max = 6.20  # Maximum chord size [m]
    chord_root = 5.38  # Chord size at the root [m]

class IIIB:
    I = 0.16

'''
Step 2: Rotor radius and rated wind speed for IIIB class
'''

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

'''
Step 4: Upscaling the blade thickness and centerline chord length
'''

# Inputs scaled to IIIB class 10 MW turbine
SF = IIIB.R / IA.R # Scaling factor Assuming linear scaling
IIIB.r_hub = SF * IA.r_hub  # Hub radius [m]
# IIIB.r = np.linspace(IIIB.r_hub, IIIB.R - 0.1, 40)  # Rotor span [m]
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
print(f"Maximum blade radius for IIIB class: {IIIB.r[-1]:.2f} m")


# Plot absolute thickness distributions
fig, ax = plt.subplots(1, 1)
ax.plot(IA.r, IA.t, label='DTU 10MW Class IA')
ax.plot(IIIB.r, IIIB.t, label='DTU 10MW Class IIIB')
ax.set_xlabel('Blade Span [m]')
ax.set_ylabel('Absolute Thickness [m]')
ax.grid(True, linestyle = ':')
ax.legend()
fig.tight_layout()
fig.savefig('A1 Aeroelastic design/Figures/absolute_thickness.svg', format='svg')

'''
Step 6: Choice of design lift and airfoils for 3 different airfoils
'''

pc_path = r'hawc_files/dtu_10mw/data/DTU_10MW_RWT_pc.dat'
tc_data = np.loadtxt(pc_path, skiprows=2, usecols=[0,1,2])
tc_241_data = tc_data[52-2:69-2, :]
tc_301_data = tc_data[158-2:175-2, :]
tc_360_data = tc_data[264-2:281-2, :]
tc_480_data = tc_data[370-2:379-2, :]

# List of tc datasets for convenience
tc_datasets = [tc_241_data, tc_301_data, tc_360_data, tc_480_data]
titles = ['tc_241_data', 'tc_301_data', 'tc_360_data', 'tc_480_data']

# Chosen design points for each airfoil
Cl_points = [1.3, 1.23, 1.37, 0.54]
Cd_points = [0.013, 0.014, 0.021, 0.032]
Cl_Cd_points = np.array(Cl_points)/np.array(Cd_points)
#Cl_Cd_points = []
#for i in range(0,len(Cl_points)):
#    Cl_Cd_points.append(Cl_points[i]/Cd_points[i])

alpha_points = [8.1, 7.5, 5.8, 1.8]
t_c_points = [24.1, 30.1, 36.0, 48.0]

# Plot for each tc dataset
for i, tc_data in enumerate(tc_datasets):
    alpha = tc_data[:, 0]  # Alpha (angle of attack)
    c_l = tc_data[:, 1]    # Lift coefficient (C_l)
    c_d = tc_data[:, 2]    # Drag coefficient (C_d)
    c_l_c_d = c_l / c_d    # Lift-to-drag ratio (C_l/C_d)

    # Create a figure for each tc_data set
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: C_l vs C_d
    axs[0].plot(c_d, c_l, 'o-', color='tab:blue')
    axs[0].set_title(f'C_l vs C_d')
    axs[0].set_xlabel('C_d')
    axs[0].set_ylabel('C_l')
    axs[0].grid(True, linestyle=':')
    axs[0].axhline(y=Cl_points[i], color='k', linestyle='--', label=f'Chosen Cl={Cl_points[i]}')
    axs[0].legend()

    # Plot 2: C_l vs Alpha
    axs[1].plot(alpha, c_l, 'o-', color='tab:orange')
    axs[1].set_title(f'C_l vs Alpha')
    axs[1].set_xlabel('Alpha [°]')
    axs[1].set_ylabel('C_l')
    axs[1].grid(True, linestyle=':')
    axs[1].axhline(y=Cl_points[i], color='k', linestyle='--', label=f'Chosen Cl={Cl_points[i]}')
    axs[1].legend()

    # Plot 3: C_l vs C_l/C_d
    axs[2].plot(c_l_c_d, c_l, 'o-', color='tab:green')
    axs[2].set_title(f'C_l vs C_l/C_d')
    axs[2].set_xlabel('C_l/C_d')
    axs[2].set_ylabel('C_l')
    axs[2].grid(True, linestyle=':')
    axs[2].axhline(y=Cl_points[i], color='k', linestyle='--', label=f'Chosen Cl={Cl_points[i]}')
    axs[2].legend()

    # Adjust layout and display the plots
    plt.tight_layout()
    plt.savefig(f'A1 Aeroelastic design/Figures/{titles[i]}.svg', format='svg')

# DTU 10MW IIIB
fig1, axs1 = plt.subplots(3, 1, figsize=(10, 8))
fig2, axs2 = plt.subplots(3, 1, figsize=(10, 8))
fig3, axs3 = plt.subplots(2, 2, figsize=(10, 8))
fig4, axs4 = plt.subplots(3, 1, figsize=(10, 8))

i = np.arange(1, 4)
colors = ['tab:blue', 'tab:orange', 'tab:green']
for i_design in i:
    # Assign the design functions for the IA class
    IA.r = np.linspace(IA.r_hub, IA.R - 0.1, 40)  # Rotor span [m]
    IA.cl_des, IA.cd_des, IA.aoa_des, IA.tc_vals, IA.cl_vals, IA.cd_vals, IA.aoa_vals = get_design_functions(i_design)

    # Design the rotor
    IA.chord, IA.tc, IA.twist, IA.cl, IA.cd, IA.aoa, IA.a, IA.CLT, IA.CLP, IA.CT, IA.CP = single_point_design(
        IA.r, IA.t, tsr, IA.R, IA.cl_des, IA.cd_des, IA.aoa_des, IA.chord_root, IA.chord_max, B)

    # Plotting design functions
    tc_plot = np.linspace(0, 100, 100)

    axs1[0].plot(tc_plot, IA.cl_des(tc_plot), color = colors[i_design-1], label=f'Design func {i_design}')
    axs1[0].plot(IA.tc_vals, IA.cl_vals, "o", color = colors[i_design-1])
    if i_design == 3:
        axs1[0].plot(t_c_points, Cl_points, 'xk', label='Design points')
    else:
        axs1[0].plot(t_c_points, Cl_points, 'xk')
    axs1[0].set_ylabel("$C_l$ [-]")
    axs1[0].set_xlim(0, 100)
    axs1[0].legend()
    axs1[0].grid(True, linestyle = ':')

    axs1[1].plot(tc_plot, IA.cl_des(tc_plot)/IA.cd_des(tc_plot), color = colors[i_design-1])
    axs1[1].plot(IA.tc_vals, np.array(IA.cl_vals)/np.array(IA.cd_vals), "o", color = colors[i_design-1])
    axs1[1].plot(t_c_points, Cl_Cd_points, 'xk')
    axs1[1].set_ylabel("$C_l/C_d$ [-]")
    axs1[1].set_xlim(0, 100)
    # axs1[1].legend()
    axs1[1].grid(True, linestyle = ':')

    axs1[2].plot(tc_plot, IA.aoa_des(tc_plot), color = colors[i_design-1])
    axs1[2].plot(IA.tc_vals, IA.aoa_vals, "o", color = colors[i_design-1])
    axs1[2].plot(t_c_points, alpha_points, 'xk')
    axs1[2].set_ylabel(r"$\alpha$ [°]")
    axs1[2].set_xlabel(r"$t/c$ [-]")
    axs1[2].set_xlim(0, 100)
    # axs1[2].legend()
    axs1[2].grid(True, linestyle = ':')

    # Plot the chord, twist and relative-thickness
    # Chord
    axs2[0].plot(IA.r, IA.chord, label=f'Design f {i_design}')
    axs2[0].set_ylabel("Chord [m]")
    axs2[0].set_xlim(0, max(IA.R, IA.R))
    axs2[0].legend()
    axs2[0].grid(True, linestyle = ':')

    # Twist
    axs2[1].plot(IA.r, IA.twist)
    axs2[1].set_ylabel(r"Twist [°]")
    axs2[1].set_xlim(0, max(IA.R, IA.R))
    # axs2[1].legend()
    axs2[1].grid(True, linestyle = ':')

    # t/c
    axs2[2].plot(IA.r, IA.tc)
    axs2[2].set_ylabel("Rel. thickness [%]")
    axs2[2].set_xlabel("Rotor span [m]")
    axs2[2].set_xlim(0, max(IA.R, IA.R))
    # axs2[2].legend()
    axs2[2].grid(True, linestyle = ':')

    # Plot r vs. t/c, aoa, cl, cd
    # t/c
    axs3[0, 0].plot(IA.r, IA.tc, label=f'Design f {i_design}')
    axs3[0, 0].set_ylabel("t/c [%]")
    axs3[0, 0].set_xlim(0, max(IA.R, IA.R))
    axs3[0, 0].legend()
    axs3[0, 0].grid(True, linestyle = ':')

    # aoa
    axs3[0, 1].plot(IA.r, IA.aoa)
    axs3[0, 1].set_ylabel(r"$\alpha$ [°]")
    axs3[0, 1].set_xlim(0, max(IA.R, IA.R))
    axs3[0, 1].yaxis.tick_right()
    axs3[0, 1].yaxis.set_label_position("right")
    # axs3[0, 1].legend()
    axs3[0, 1].grid(True, linestyle = ':')

    # cl
    axs3[1, 0].plot(IA.r, IA.cl)
    axs3[1, 0].set_ylabel("$C_l$ [-]")
    axs3[1, 0].set_xlabel("Span [m]")
    axs3[1, 0].set_xlim(0, max(IA.R, IA.R))
    # axs3[1, 0].legend()
    axs3[1, 0].grid(True, linestyle = ':')

    # cd
    axs3[1, 1].plot(IA.r, IA.cd)
    axs3[1, 1].set_ylabel("$C_d$ [-]")
    axs3[1, 1].set_xlabel("Span [m]")
    axs3[1, 1].set_xlim(0, max(IA.R, IA.R))
    axs3[1, 1].yaxis.tick_right()
    axs3[1, 1].yaxis.set_label_position("right")
    # axs3[1, 1].legend()
    axs3[1, 1].grid(True, linestyle = ':')

    # Plot r vs. CLT, CLP, a
    # Local-Thrust-Coefficient
    axs4[0].plot(IA.r, IA.CLT, label=f'Design f {i_design}')
    axs4[0].axhline(y=8 / 9, ls="--", color="k")
    axs4[0].set_ylabel("Local thrust ($C_{LT}$) [-]")
    axs4[0].set_ylim(0, 1.0)
    axs4[0].set_xlim(0, max(IA.R, IA.R))
    axs4[0].legend()
    axs4[0].grid(True, linestyle = ':')

    # Local-Power-Coefficient
    axs4[1].plot(IA.r, IA.CLP)
    axs4[1].axhline(y=16 / 27, ls="--", color="k")
    axs4[1].set_ylabel("Local Power ($C_{LP}$) [-]")
    axs4[1].set_xlim(0, max(IA.R, IA.R))
    axs4[1].set_ylim(-0.4, 0.6)
    # axs4[1].legend()
    axs4[1].grid(True, linestyle = ':')

    # AxIAl Induction
    axs4[2].plot(IA.r, IA.a)
    axs4[2].axhline(y=1 / 3, ls="--", color="k")
    axs4[2].set_ylabel("AxIAl induction ($a$) [-]")
    axs4[2].set_xlabel("Rotor span [m]")
    axs4[2].set_xlim(0, max(IA.R, IA.R))
    # axs4[2].legend()
    axs4[2].grid(True, linestyle = ':')

    fig4.suptitle(f"$C_T$={IA.CT:1.3f}, $C_P$={IA.CP:1.3f}")
    
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig1.savefig('A1 Aeroelastic design/Figures/design_functions_clcd.png', format='png')
fig1.savefig('A1 Aeroelastic design/Figures/design_functions_clcd.svg', format='svg')
fig2.savefig('A1 Aeroelastic design/Figures/chord_twist_thickness.svg', format='svg')
fig3.savefig('A1 Aeroelastic design/Figures/aoa_cl_cd.svg', format='svg')
fig4.savefig('A1 Aeroelastic design/Figures/CLT_CLP_a.svg', format='svg')

# Design function 3 is chosen
i_design = 3

'''
Step 7: Find the design with tip-speed-ratio that maximizes CP
'''

tsr_range = np.arange(6, 12, 1)
IA.CP_store = np.zeros(len(tsr_range))
IA.cl_des, IA.cd_des, IA.aoa_des, IA.tc_vals, IA.cl_vals, IA.cd_vals, IA.aoa_vals = get_design_functions(3)

for i, tsr in enumerate(tsr_range):
    IA.chord, IA.tc, IA.twist, IA.cl, IA.cd, IA.aoa, IA.a, IA.CLT, IA.CLP, IA.CT, IA.CP = single_point_design(
        IA.r, IA.t, tsr, IA.R, IA.cl_des, IA.cd_des, IA.aoa_des, IA.chord_root, IA.chord_max, B)
    IA.CP_store[i] = IA.CP

    print(f"tsr: {tsr}, CP: {IA.CP:1.3f}")

# Plot the power coefficient
plt.figure()
plt.plot(tsr_range, IA.CP_store, linestyle = '-', marker = 'o')
plt.axhline(y=max(IA.CP_store), color='grey', linestyle='--')
plt.axvline(x=tsr_range[np.argmax(IA.CP_store)], color='grey', linestyle='--')
plt.xlabel("Tip-Speed-Ratio [-]")
plt.ylabel("Power Coefficient [-]")
plt.grid(True, linestyle = ':')
plt.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures/CP_vs_TSR.svg', format='svg')

# TSR 7 is chosen
tsr = 7

'''
Step 8: Present your final chord, twist, and relative thickness distributions and compare to the original DTU 10MW rotor. 
'''

# Final designs for chosen TSR = 7 and design function 3
# IA class
IA.cl_des, IA.cd_des, IA.aoa_des, IA.tc_vals, IA.cl_vals, IA.cd_vals, IA.aoa_vals = get_design_functions(i_design)

IA.chord, IA.tc, IA.twist, IA.cl, IA.cd, IA.aoa, IA.a, IA.CLT, IA.CLP, IA.CT, IA.CP = single_point_design(
    IA.r, IA.t, tsr, IA.R, IA.cl_des, IA.cd_des, IA.aoa_des, IA.chord_root, IA.chord_max, B)

# IIB class
IIIB.r = np.linspace(IIIB.r_hub, IIIB.R - 0.1, 40)
IIIB.cl_des, IIIB.cd_des, IIIB.aoa_des, IIIB.tc_vals, IIIB.cl_vals, IIIB.cd_vals, IIIB.aoa_vals = get_design_functions(i_design)

IIIB.chord, IIIB.tc, IIIB.twist, IIIB.cl, IIIB.cd, IIIB.aoa, IIIB.a, IIIB.CLT, IIIB.CLP, IIIB.CT, IIIB.CP = single_point_design(
        IIIB.r, IIIB.t, tsr, IIIB.R, IIIB.cl_des, IIIB.cd_des, IIIB.aoa_des, IIIB.chord_root, IIIB.chord_max, B)

# Plot the chord, twist and relative-thickness of the IA and IIIB classes
fig, axs = plt.subplots(3, 1, figsize=(10, 8))

# Chord
axs[0].plot(IA.r, IA.chord, label='DTU 10MW Class IA')
axs[0].plot(IIIB.r, IIIB.chord, label='DTU 10MW Class IIIB')
axs[0].set_ylabel("Chord [m]")
axs[0].legend()
axs[0].grid(True, linestyle = ':')

# Twist
axs[1].plot(IA.r, IA.twist, label='DTU 10MW Class IA')
axs[1].plot(IIIB.r, IIIB.twist, label='DTU 10MW Class IIIB')
axs[1].set_ylabel(r"Twist [°]")
# axs[1].legend()
axs[1].grid(True, linestyle = ':')

# Absolute thickness
axs[2].plot(IA.r, IA.t, label='DTU 10MW Class IA')
axs[2].plot(IIIB.r, IIIB.t, label='DTU 10MW Class IIIB')
axs[2].set_ylabel("Absolute Thickness [m]")
axs[2].set_xlabel("Blade Span [m]")
# axs[2].legend()
axs[2].grid(True, linestyle = ':')

plt.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures/final_chord_twist_thickness.svg', format='svg')



# plt.show()


# %%
# --------Calculations for PART 3-------

#interpolate twist values from our design to the htcfile
print(IA.twist)
print(np.size(IA.twist))
old_twist = [-14.5, -14.5, -14.4851, -14.461, -14.3388, -14.0201, -13.3904, 
             -12.9371, -11.9445, -9.98243, -8.45147, -7.46417, -6.72916, 
             -6.08842, -5.49322, -4.39222, -3.09315, -1.75629, -0.50065, 
             0.60196, 1.5556, 2.51935, 2.7295, 2.93201, 3.11874, 3.28847, 3.42796]
new_twist_interpolated = np.interp(np.arange(0,27), np.arange(0,len(IA.twist)), IA.twist)

print(new_twist_interpolated)

TSR = 7.5  # Tip speed ratio
R = 89.5  # Radius in meters
v = 11.13  # Wind speed in meters per second
# Calculate angular velocity (omega) in rad/s
omega = TSR * v / R
# Convert angular velocity to RPM
RPM = (omega * 60) / (2 * np.pi)

print(f"Max RPM: {RPM}")
print(f"Max generator speed given a gear ratio of 50: {RPM*50}")

# save data from design part to file
# Save the arrays to a file in the current directory
np.savez('3IIIB_design_data.npz', r=IIIB.r, cl=IIIB.cl(IIIB.r), cd=IIIB.cd(IIIB.r), aoa=IIIB.aoa(IIIB.r))