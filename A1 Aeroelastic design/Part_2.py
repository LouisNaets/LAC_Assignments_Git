# %% Import modules
import matplotlib.pyplot as plt
import numpy as np
from aero_design_functions import get_design_functions, single_point_design
from sympy import symbols, Eq, solve
from matplotlib import cm
from lacbox.io import load_pwr, load_ind, load_inds, load_st, save_st
from structural_scaling_example import scale_ST_data

"""
Assignment 1: Aeroelastic Design part 2
DTU 10MW Class IA -> IIIB
"""

# Global constants
rho = 1.225 # Air Density [kg/m^3]
P_rated = 10e6 # Rated Power [W]
tsr = 9.0 # Tip-Speed-Ratio [-]
B = 3  # Number of blades [#]
V_0 = 8 # Incoming wind speed [m/s]

# IEC 61400-1 Classes
class IA:
    I = 0.18 # Turbulence Intensity[-]
    V_rated = 11.4 # m/s 
    r_hub = 2.8  # Hub radius [m]
    R = 89.15 # Rotor diameter (excl hub) [m]
    chord_max = 6.20  # Maximum chord size [m]
    chord_root = 5.38  # Chord size at the root [m]

class IIIB:
    I = 0.16

'''
Step 2: Rotor radius and rated wind speed for IIIB class
'''

plt.rcParams.update({'axes.labelsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12})  # Affects both X and Y labels

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
IIIB.chord_root = IA.chord_root  # !! shouldn't be scaled !!

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
#IIIB.tc = SF * IA.tc
IIIB.t = IA.t # Absolute Thickness [m] !!! should only be scaled in the x-direction, so the change in radius should just be adjusted instead of using SF here !!!

print(f"Chord size at the root for IIIB class: {IIIB.chord_root:.2f} m")
print(f"Maximum chord size for IIIB class: {IIIB.chord_max:.2f} m")

# Plot absolute thickness distributions
fig, ax = plt.subplots(1, 1)
ax.plot(IA.r, IA.t, label='DTU 10MW rotor')
ax.plot(IIIB.r, IIIB.t, label='New IIIB design')
ax.set_xlabel('Blade span [m]')
ax.set_ylabel('Absolute thickness [m]')
ax.grid(True, linestyle = ':')
ax.legend()
fig.tight_layout()
fig.savefig('A1 Aeroelastic design/Figures/absolute_thickness.svg', format='svg')
fig.savefig('A1 Aeroelastic design/Figures/absolute_thickness.png', format='png')

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
    #axs[0].set_title(f'C_l vs C_d')
    axs[0].set_xlabel(r'$C_d$ [-]', fontsize=16)
    axs[0].set_ylabel(r'$C_l$ [-]', fontsize=16)
    axs[0].grid(True, linestyle=':')
    axs[0].axhline(y=Cl_points[i], color='k', linestyle='--', label=f'Chosen $C_l$={Cl_points[i]}')
    axs[0].legend(fontsize=14)

    # Plot 2: C_l vs Alpha
    axs[1].plot(alpha, c_l, 'o-', color='tab:orange')
    #axs[1].set_title(f'C_l vs Alpha')
    axs[1].set_xlabel(r'$\alpha$ [°]', fontsize=16)
    axs[1].set_ylabel(r'$C_l$ [-]', fontsize=16)
    axs[1].grid(True, linestyle=':')
    axs[1].axhline(y=Cl_points[i], color='k', linestyle='--', label=f'Chosen $C_l$={Cl_points[i]}')
    axs[1].legend(fontsize=14)

    # Plot 3: C_l vs C_l/C_d
    axs[2].plot(c_l_c_d, c_l, 'o-', color='tab:green')
    #axs[2].set_title(f'C_l vs C_l/C_d')
    axs[2].set_xlabel(r'$C_l$/$C_d$ [-]', fontsize=16)
    axs[2].set_ylabel(r'$C_l$ [-]', fontsize=16)
    axs[2].grid(True, linestyle=':')
    axs[2].axhline(y=Cl_points[i], color='k', linestyle='--', label=f'Chosen $C_l$={Cl_points[i]}')
    axs[2].legend(fontsize=14)

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
    axs1[2].set_xlabel(r"Relative thickness [%]")
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
    axs2[2].set_ylabel("Relative thickness [%]")
    axs2[2].set_xlabel("Rotor span [m]")
    axs2[2].set_xlim(0, max(IA.R, IA.R))
    # axs2[2].legend()
    axs2[2].grid(True, linestyle = ':')

    # Plot r vs. t/c, aoa, cl, cd
    # t/c
    axs3[0, 0].plot(IA.r, IA.tc, label=f'Design f {i_design}')
    axs3[0, 0].set_ylabel("Relative thickness [%]")
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

tsr_range = np.arange(6, 12, 0.5)
IA.CP_store = np.zeros(len(tsr_range))
IA.cl_des, IA.cd_des, IA.aoa_des, IA.tc_vals, IA.cl_vals, IA.cd_vals, IA.aoa_vals = get_design_functions(i_design)

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
plt.xlabel("Tip-speed ratio (TSR) [-]")
plt.ylabel(r"$C_p$ [-]")
plt.grid(True, linestyle = ':')
plt.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures/CP_vs_TSR.svg', format='svg')
plt.savefig('A1 Aeroelastic design/Figures/CP_vs_TSR.png', format='png')

# TSR 7.5 is chosen
tsr = 7.5

'''
Step 8: Present your final chord, twist, and relative thickness distributions and compare to the original DTU 10MW rotor. 
'''

# Final designs for chosen TSR = 7.5 and design function 3

# !! This wasn't actually required we just needed the original values for the DTU 10MW reference turbine
# IA class 
#IA.cl_des, IA.cd_des, IA.aoa_des, IA.tc_vals, IA.cl_vals, IA.cd_vals, IA.aoa_vals = get_design_functions(i_design)

#IA.chord, IA.tc, IA.twist, IA.cl, IA.cd, IA.aoa, IA.a, IA.CLT, IA.CLP, IA.CT, IA.CP = single_point_design(
#    IA.r, IA.t, tsr, IA.R, IA.cl_des, IA.cd_des, IA.aoa_des, IA.chord_root, IA.chord_max, B)


#From DTU 10MW master htc file
IA.twist_r = [4.44089E-16,3.00000E+00,6.00000E+00,7.00004E+00,8.70051E+00,1.04020E+01,1.22046E+01,1.32065E+01,
            1.50100E+01,1.82151E+01,2.14178E+01,2.46189E+01,2.78193E+01,3.10194E+01,3.42197E+01,4.02204E+01,
            4.66217E+01,5.30232E+01,5.94245E+01,6.58255E+01,7.22261E+01,7.90266E+01,8.05267E+01,8.20271E+01,
            8.35274E+01,8.50277E+01,8.63655E+01]
IA.twist = -1*np.array([-1.45000E+01,-1.45000E+01,-1.44851E+01,-1.44610E+01,-1.43388E+01,-1.40201E+01,-1.33904E+01,
              -1.29371E+01,-1.19445E+01,-9.98243E+00,-8.45147E+00,-7.46417E+00,-6.72916E+00,-6.08842E+00,
              -5.49322E+00,-4.39222E+00,-3.09315E+00,-1.75629E+00,-5.00650E-01,6.01964E-01,1.55560E+00,
              2.51935E+00,2.72950E+00,2.93201E+00,3.11874E+00,3.28847E+00,3.42796E+00])

# IIIB class
IIIB.r = np.linspace(IIIB.r_hub, IIIB.R - 0.1, 40)
IIIB.cl_des, IIIB.cd_des, IIIB.aoa_des, IIIB.tc_vals, IIIB.cl_vals, IIIB.cd_vals, IIIB.aoa_vals = get_design_functions(i_design)

IIIB.chord, IIIB.tc, IIIB.twist, IIIB.cl, IIIB.cd, IIIB.aoa, IIIB.a, IIIB.CLT, IIIB.CLP, IIIB.CT, IIIB.CP = single_point_design(
        IIIB.r, IIIB.t, tsr, IIIB.R, IIIB.cl_des, IIIB.cd_des, IIIB.aoa_des, IIIB.chord_root, IIIB.chord_max, B)

# Plot the chord, twist and relative-thickness of the DTU 10MW reference turbine and new IIIB design
# Determine chord and thickness for IA class from the _ae file
file_path = r'hawc_files\dtu_10mw\data\DTU_10MW_RWT_ae.dat' 
data = np.loadtxt(file_path, skiprows=2, usecols=(0, 1, 2))
IA.r = data[:, 0]  # Blade span [m] - please note: the ae file starts after the hub
IA.chord = data[:, 1]  # Chord length [m]
IA.tc = data[:, 2]  # Thickness-to-chord ratio [-]
IA.t = IA.tc/100 * IA.chord  # Absolute Thickness [m]

fig, axs = plt.subplots(3, 1, figsize=(10, 8))

# Chord
axs[0].plot(IA.r, IA.chord, label='DTU 10MW rotor')
axs[0].plot(IIIB.r-IIIB.r_hub, IIIB.chord, label='New IIIB design')
axs[0].set_ylabel("Chord [m]")
axs[0].legend()
axs[0].grid(True, linestyle = ':')

# Twist
axs[1].plot(IA.twist_r, IA.twist, label='DTU 10MW rotor')
axs[1].plot(IIIB.r-IIIB.r_hub, IIIB.twist, label='New IIIB design')
axs[1].set_ylabel(r"Twist [°]")
# axs[1].legend()
axs[1].grid(True, linestyle = ':')

# Relative thickness
axs[2].plot(IA.r, IA.tc, label='DTU 10MW rotor')
axs[2].plot(IIIB.r-IIIB.r_hub, IIIB.tc, label='New IIIB design')
axs[2].set_ylabel("Relative thickness [%]")
axs[2].set_xlabel("Blade span [m]")
# axs[2].legend()
axs[2].grid(True, linestyle = ':')

plt.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures/final_chord_twist_thickness.png', format='png')
plt.savefig('A1 Aeroelastic design/Figures/final_chord_twist_thickness.svg', format='svg')
#plt.show()


# %%

#Part 3

IA.htc_main_z = [4.44089E-16, 3.00000E+00, 6.00000E+00, 7.00004E+00, 8.70051E+00,
1.04020E+01, 1.22046E+01, 1.32065E+01, 1.50100E+01, 1.82151E+01, 2.14178E+01,
2.46189E+01, 2.78193E+01, 3.10194E+01, 3.42197E+01, 4.02204E+01, 4.66217E+01,
5.30232E+01, 5.94245E+01, 6.58255E+01, 7.22261E+01, 7.90266E+01, 8.05267E+01,
8.20271E+01, 8.35274E+01, 8.50277E+01, 8.63655E+01]

IIIB.htc_main_x = [0.00000E+00, -2.06477E-05, -7.28810E-03, -1.89235E-02, -5.41282E-02,
-1.26633E-01, -2.25666E-01, -2.88563E-01, -3.99194E-01, -5.76634E-01, -7.07136E-01,
-7.91081E-01, -8.37195E-01,	-8.53948E-01, -8.49367E-01,	-7.93920E-01, -7.16284E-01,
-6.34358E-01, -5.53179E-01,	-4.75422E-01, -4.03180E-01,	-3.30085E-01, -3.10140E-01,
-2.86719E-01, -2.55823E-01,	-2.07891E-01, -8.98940E-02]

IIIB.htc_main_y = [7.00600E-05, -1.22119E-02, -2.49251E-02, -2.73351E-02, -2.82163E-02,
-2.13210E-02, -1.28378E-02, -7.70659E-03, -4.88317E-03, -1.80296E-02, -5.01772E-02,
-9.41228E-02, -1.48880E-01,	-2.14514E-01, -2.90618E-01,	-4.62574E-01, -6.88437E-01,
-9.60017E-01, -1.28424E+00,	-1.66402E+00, -2.10743E+00,	-2.65630E+00, -2.78882E+00,
-2.92517E+00, -3.06577E+00, -3.20952E+00, -3.33685E+00]

#These are used to udpate the htc main file:

IIIB.htc_main_z = np.array(IA.htc_main_z)*SF
IIIB.htc_main_twist = np.interp(IIIB.htc_main_z, IIIB.r, IIIB.twist)

#!!! haven't updated past this point after presentation !!!

print_extra = 1
if print_extra:
    print('Input for the htc main file:')
    for i in range(0, 27):
        print('sec\t' + str(i+1) + '\t' + str(round(IIIB.htc_main_x[i], 7)) + '\t' + str(round(IIIB.htc_main_y[i], 7)) + '\t' + str(round(IIIB.htc_main_z[i], 7)) + '\t' + str(round(-IIIB.htc_main_twist[i], 7)) + '\t;')

    print('Input for the ae.dat file:')
    for i in range(0,len(IIIB.r)):
        #This if-statement is required because from the ae.dat file we need a point a 0m, therefor this assumption is used.
        #if i == 0:
        #    print('0\t' + str(round(IIIB.chord[i], 7)) + '\t100\t1\t;')
        #else:
        #    print(str(round(IIIB.r[i-1], 7)) + '\t' + str(round(IIIB.chord[i-1], 7)) + '\t' + str(round(IIIB.t[i-1]/IIIB.chord[i-1]*100, 7)) + '\t1\t;')
        
        #New version - remember to remove r_hub for ae file
        print(str(round(IIIB.r[i]-IIIB.r_hub, 7)) + '\t' + str(round(IIIB.chord[i], 7)) + '\t' + str(max(round(IIIB.t[i]/IIIB.chord[i]*100, 7), 24.1)) + '\t1\t;')


print('Omegas for the opt files')
for tsr_opt in np.arange(6, 10.1, 0.5):
    omega = tsr_opt*V_0/IIIB.R * (60/(2*np.pi))
    print('TSR: ' + str(tsr_opt) + '-> omega: ' + str(round(omega,7)) + ' RPM')


#Part 3 plots 1 and 2:

# Function to plot with color map and dotted lines
def plot_with_colormap(ax, x, y, label, cmap, idx, total_lines):
    colors = cm.viridis(np.linspace(0, 1, total_lines))  # Gradual colormap
    ax.plot(x, y, label=label, color=colors[idx], linestyle='-') 

#plot_with_colormap(axes[0, 0], r_R, c_R, label=f'TSR = {TSR}', cmap=cm.viridis, idx=idx, total_lines=len(TSR_range))
#plot_with_colormap(axes[0, 1], r_R, np.rad2deg(theta_R), label=f'TSR = {TSR}', cmap=cm.viridis, idx=idx, total_lines=len(TSR_range))




# Path for the file
ind_path = "./hawc_files/our_design/res_hawc2s/group7_3B_design_hawc2s_1wsp_u8000.ind"
# Load the data
ind_data = load_ind(ind_path)

IIIB.relative_t = np.array(IIIB.t/IIIB.chord)*100
IIIB.relative_t_ind = np.interp(ind_data["s_m"], IIIB.r, IIIB.relative_t)

tc_plot_r = np.interp(tc_plot, IIIB.relative_t_ind, ind_data["s_m"]) #this is not correct yet

#print(tc_plot_r)

#print(tc_plot)
#print(IIIB.relative_t)
#print(ind_data["s_m"])
#print(IIIB.relative_t_ind)

'''FIGURE 3.1'''

plt.rcParams.update({'axes.labelsize': 18, 'xtick.labelsize': 18, 'ytick.labelsize': 18, 'legend.fontsize': 15})

fig5, axes5 = plt.subplots(3, 2, figsize=(18, 12), dpi=500)

axes5[0,0].plot(IIIB.tc, IIIB.cl, color = colors[0], label='Design $C_l$')
axes5[0,0].plot(IIIB.relative_t_ind, ind_data["Cl"], color = colors[1], label='HAWC2S $C_l$')
axes5[0,0].set_ylabel("$C_l$ [-]")
axes5[0,0].set_xlim(0, 100)
axes5[0,0].legend()
axes5[0,0].grid(True, linestyle = ':')

axes5[1,0].plot(IIIB.tc, IIIB.cl/IIIB.cd, color = colors[0], label='Design $C_l/C_d$')
axes5[1,0].plot(IIIB.relative_t_ind, ind_data["Cl"]/ind_data["Cd"], color = colors[1], label='HAWC2S $C_l/C_d$')
axes5[1,0].set_ylabel("$C_l/C_d$ [-]")
axes5[1,0].set_xlim(0, 100)
axes5[1,0].legend()
axes5[1,0].grid(True, linestyle = ':')

axes5[2,0].plot(IIIB.tc, IIIB.aoa, color = colors[0], label=r'Design $\alpha$')
axes5[2,0].plot(IIIB.relative_t_ind, np.rad2deg(ind_data["aoa_rad"]), color = colors[1], label=r'HAWC2S $\alpha$')
axes5[2,0].set_ylabel(r"$\alpha$ [°]")
axes5[2,0].set_xlabel(r"Relative thickness [%]")
axes5[2,0].set_xlim(0, 100)
axes5[2,0].legend()
axes5[2,0].grid(True, linestyle = ':')

axes5[0,1].plot(IIIB.r, IIIB.cl, color = colors[0], label='Design $C_l$')
axes5[0,1].plot(ind_data["s_m"], ind_data["Cl"], color = colors[1], label='HAWC2S $C_l$')
axes5[0,1].set_ylabel("$C_l$ [-]")
axes5[0,1].set_xlim(0, 100)
axes5[0,1].legend()
axes5[0,1].grid(True, linestyle = ':')

axes5[1,1].plot(IIIB.r, IIIB.cl/IIIB.cd, color = colors[0], label='Design $C_l/C_d$')
axes5[1,1].plot(ind_data["s_m"], ind_data["Cl"]/ind_data["Cd"], color = colors[1], label='HAWC2S $C_l/C_d$')
axes5[1,1].set_ylabel("$C_l/C_d$ [-]")
axes5[1,1].set_xlim(0, 100)
axes5[1,1].legend()
axes5[1,1].grid(True, linestyle = ':')

axes5[2,1].plot(IIIB.r, IIIB.aoa, color = colors[0], label=r'Design $\alpha$')
axes5[2,1].plot(ind_data["s_m"], np.rad2deg(ind_data["aoa_rad"]), color = colors[1], label=r'HAWC2S $\alpha$')
axes5[2,1].set_ylabel(r"$\alpha$ [°]")
axes5[2,1].set_xlabel(r"Curvelinear radius [m]")
axes5[2,1].set_xlim(0, 100)
axes5[2,1].legend()
axes5[2,1].grid(True, linestyle = ':')
fig5.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part3/3.1.svg', format='svg')

# Path for the file
inds_path = "./hawc_files/our_design/res_hawc2s/group7_3B_design_hawc2s_multitsr_u8003.ind"
# Load the data
inds_data = load_ind(inds_path)

# Path for the file
pwr_path = "./hawc_files/our_design/res_hawc2s/group7_3B_design_hawc2s_multitsr.pwr"
# Load the data
pwr_data = load_pwr(pwr_path)



'''FIGURE 3.2'''

# Set up a 3x2 grid layout
fig6, axes6 = plt.subplots(3, 2, figsize=(18, 12), dpi=500)

axes6[0,0].plot(inds_data["s_m"], inds_data["Cl"], color = colors[0])
axes6[0,0].set_ylabel("$C_l$ [-]")
axes6[0,0].grid(True, linestyle = ':')

axes6[0,1].plot(inds_data["s_m"], inds_data["Cl"]/inds_data["Cd"], color = colors[0])
axes6[0,1].set_ylabel("$C_l$/$C_d$ [-]")
axes6[0,1].grid(True, linestyle = ':')

axes6[1,0].plot(inds_data["s_m"], np.rad2deg(inds_data["aoa_rad"]), color = colors[0])
axes6[1,0].set_ylabel(r"$\alpha$ [°]")
axes6[1,0].grid(True, linestyle = ':')

axes6[1,1].plot(inds_data["s_m"], inds_data["a"], color = colors[0])
axes6[1,1].axhline(y=1/3, color='black', linestyle='--', label='Optimum $a$')
axes6[1,1].set_ylabel(r"Axial induction [-]")
axes6[1,1].grid(True, linestyle = ':')
axes6[1,1].legend()

axes6[2,0].plot(inds_data["s_m"], inds_data["CP"], color = colors[0])
axes6[2,0].axhline(y=16/27, color='black', linestyle='--', label='Betz limit')
axes6[2,0].set_ylabel(r"Local $C_p$ [-]")
axes6[2,0].set_xlabel(r"Curvelinear radius [m]")
axes6[2,0].grid(True, linestyle = ':')
axes6[2,0].legend()

axes6[2,1].plot(inds_data["s_m"], inds_data["CT"], color = colors[0])
axes6[2,1].set_ylabel(r"Local $C_T$ [-]")
axes6[2,1].set_xlabel(r"Curvelinear radius [m]")
axes6[2,1].grid(True, linestyle = ':')

# Save the figure
fig6.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part3/3.2.svg', format='svg')

plt.rcParams.update({'axes.labelsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12})

'''FIGURE 3.3'''

# Set up a 1x2 grid layout
tsrs = np.arange(6, 10.1, 0.5)

fig7, axes7 = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

axes7[0].plot(tsrs, pwr_data["Cp"], color = colors[0])
axes7[0].set_ylabel(r"$C_p$ [-]")
axes7[0].set_xlabel(r"Tip-speed ratio [-]")
axes7[0].grid(True, linestyle = ':')

axes7[1].plot(tsrs, pwr_data["Ct"], color = colors[0])
axes7[1].set_ylabel(r"$C_T$ [-]")
axes7[1].set_xlabel(r"Tip-speed ratio [-]")
axes7[1].grid(True, linestyle = ':')

# Save the figure
fig7.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part3/3.3.svg', format='svg')

'''FIGURE 3.4'''

opt_path = './hawc_files/our_design/data/group7_3B_design_rigid.opt'
rigid_opt_data = np.loadtxt(opt_path, skiprows=1)

class r_opt:
    V_o = rigid_opt_data[:, 0]
    pitch = rigid_opt_data[:, 1]
    omega = rigid_opt_data[:, 2]
    P = rigid_opt_data[:, 3]
    T = rigid_opt_data[:, 4]
    C_p = []
    C_T = []

# Set up a 1x2 grid layout
fig8, axes8 = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

axes8[0].plot(r_opt.V_o, r_opt.omega, color = colors[0])
axes8[0].axvline(x=IIIB.V_rated, color='black', linestyle='--', label='Rated wind speed')
axes8[0].set_ylabel(r"$\omega$ [RPM]")
axes8[0].set_xlabel(r"Wind speed [m/s]")
axes8[0].grid(True, linestyle = ':')
axes8[0].legend()

axes8[1].plot(r_opt.V_o, r_opt.pitch, color = colors[0])
axes8[1].axvline(x=IIIB.V_rated, color='black', linestyle='--', label='Rated wind speed')
axes8[1].set_ylabel(r"$\theta_p$ [°]")
axes8[1].set_xlabel(r"Wind speed [m/s]")
axes8[1].grid(True, linestyle = ':')
axes8[1].legend()

# Save the figure
fig8.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part3/3.4.svg', format='svg')



'''FIGURE 3.5'''

rho = 1.225 #kg/m3
area = IIIB.R**2*np.pi
for i in range(0,len(r_opt.V_o)):
    r_opt.C_p.append(r_opt.P[i]*1000/(0.5*rho*area*r_opt.V_o[i]**3))
    r_opt.C_T.append(r_opt.T[i]*1000/(0.5*rho*area*r_opt.V_o[i]**2))


# Set up a 1x2 grid layout
fig9, axes9 = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

# Combined plot for Power and Thrust
ax1 = axes9[0]

# Plot Power
power_line = ax1.plot(r_opt.V_o, r_opt.P, color=colors[0], label='Aerodynamic power')[0]
#ax1.axvline(x=IIIB.V_rated, color='black', linestyle='--', label='Rated wind speed')

# Create a twin axis for Thrust
ax2 = ax1.twinx()
thrust_line = ax2.plot(r_opt.V_o, r_opt.T, color=colors[1], label='Thrust')[0]

# Set labels
ax1.set_ylabel(r"Power [kW]", color=colors[0])
ax1.set_xlabel(r"Wind speed [m/s]")
ax1.grid(True, linestyle=':')
ax1.tick_params(axis='y', labelcolor=colors[0])
ax2.set_ylabel(r"Thrust [kN]", color=colors[1])
ax2.tick_params(axis='y', labelcolor=colors[1])

# Combine legends for Power and Thrust
ax1.legend(handles=[power_line, thrust_line])

# Plot for C_p and C_T
ax3 = axes9[1]

# Plot C_p
cp_line = ax3.plot(r_opt.V_o, r_opt.C_p, color=colors[0], label=r'$C_p$')[0]
#ax3.axvline(x=IIIB.V_rated, color='black', linestyle='--', label='Rated wind speed')

# Create a twin axis for C_T
ax4 = ax3.twinx()
ct_line = ax4.plot(r_opt.V_o, r_opt.C_T, color=colors[1], label=r'$C_T$')[0]

# Set labels
ax3.set_ylabel(r"Power coefficient [-]", color=colors[0])
ax3.set_xlabel(r"Wind speed [m/s]")
ax3.grid(True, linestyle=':')
ax3.tick_params(axis='y', labelcolor=colors[0])
ax4.set_ylabel(r"Thrust coefficient [-]", color=colors[1])
ax4.tick_params(axis='y', labelcolor=colors[1])

ax3.set_ylim(0, 0.85)
ax4.set_ylim(0, 0.85)

# Combine legends for C_p and C_T
ax3.legend(handles=[cp_line, ct_line])

# Adjust layout and save the figure
fig9.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part3/3.5.svg', format='svg')

'''
PART 4
'''

path_st_file_DTU10MW = "./hawc_files/our_design/data/DTU_10MW_RWT_Blade_st.dat"
st_data_DTU10MW = load_st(path_st_file_DTU10MW, 0, 0)  # Baseline data

new_ST_data = scale_ST_data(st_data_DTU10MW, SF)

save_st("./hawc_files/our_design/data/group7_3B_design_Blade_st.dat", new_ST_data)

opt_path = './hawc_files/our_design/data/group7_3B_design_flex.opt'
flex_opt_data = np.loadtxt(opt_path, skiprows=1)

class f_opt:
    V_o = flex_opt_data[:, 0]
    pitch = flex_opt_data[:, 1]
    omega = flex_opt_data[:, 2]
    P = flex_opt_data[:, 3]
    T = flex_opt_data[:, 4]
    C_p = []
    C_T = []

for i in range(0,len(f_opt.V_o)):
    f_opt.C_p.append(f_opt.P[i]*1000/(0.5*rho*area*f_opt.V_o[i]**3))
    f_opt.C_T.append(f_opt.T[i]*1000/(0.5*rho*area*f_opt.V_o[i]**2))

Blade_flex_st_path = './hawc_files/our_design/data/group7_3B_design_Blade_st.dat'
Blade_flex_st_data = np.loadtxt(Blade_flex_st_path, skiprows=3)

class blade_st_dat:
    def __init__(self, Blade_flex_st_data):
        self.r = Blade_flex_st_data[:, 0]  # radius
        self.distMass = Blade_flex_st_data[:, 1]  # distributed mass [kg/m]
        self.I_x = Blade_flex_st_data[:, 10]  # area moment of inertia about x-axis
        self.I_y = Blade_flex_st_data[:, 11]  # area moment of inertia about y-axis
        self.crossArea = Blade_flex_st_data[:, 15]  # cross-sectional area [m^2]
        
        # Calculate mass moment of inertia
        self.I_x_mass = self.calculate_mass_moments(self.I_x)
        self.I_y_mass = self.calculate_mass_moments(self.I_y)
    
    def calculate_mass_moments(self, I_area):
        # Assuming uniform segment length (change this if necessary)
        segment_length = np.diff(self.r)
        distMassAve = (self.distMass[:-1] + self.distMass[1:]) / 2
        masses = distMassAve * segment_length  # Total mass in each segment
        
        # Calculate mass moments of inertia
        I_area_ave = (I_area[:-1] + I_area[1:]) / 2
        crossArea_ave = (self.crossArea[:-1] + self.crossArea[1:]) / 2

        I_mass = I_area_ave * masses  / crossArea_ave# Convert area moments to mass moments
        return I_mass

flex_blade_st_dat = blade_st_dat(Blade_flex_st_data)

'''FIGURE 4.1'''

# Set up a 1x2 grid layout
fig10, axes10 = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

axes10[0].plot(flex_blade_st_dat.r, flex_blade_st_dat.distMass, color = colors[0])
axes10[0].set_ylabel(r"Distributed mass [kg/m]")
axes10[0].set_xlabel(r"Curvelinear radius [m]")
axes10[0].grid(True, linestyle = ':')

axes10[1].plot(flex_blade_st_dat.r[1:], flex_blade_st_dat.I_x_mass, color = colors[0], label = r'$I_x$')
axes10[1].plot(flex_blade_st_dat.r[1:], flex_blade_st_dat.I_y_mass, color = colors[1], label = r'$I_y$')
axes10[1].set_ylabel(r"Mass moment of inertia [kg$\cdot$m$^2$]")
axes10[1].set_xlabel(r"Curvelinear radius [m]")
axes10[1].legend()
axes10[1].grid(True, linestyle = ':')

# Save the figure
fig10.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part4/4.1.png', format='png')



'''FIGURE 4.2'''

# Set up a 1x2 grid layout
fig11, axes11 = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

axes11[0].plot(r_opt.V_o, r_opt.C_p, color = colors[0], label='Rigid')
axes11[0].plot(f_opt.V_o, f_opt.C_p, color = colors[1], linestyle='--', label='Flexible')
axes11[0].set_ylabel(r"Power coefficient [-]")
axes11[0].set_xlabel(r"Wind speed [m/s]")
axes11[0].legend()
axes11[0].grid(True, linestyle = ':')

axes11[1].plot(r_opt.V_o, r_opt.C_T, color = colors[0], label='Rigid')
axes11[1].plot(f_opt.V_o, f_opt.C_T, color = colors[1], linestyle='--', label='Flexible')
axes11[1].set_ylabel(r"Thrust coefficient [-]")
axes11[1].set_xlabel(r"Wind speed [m/s]")
axes11[1].legend()
axes11[1].grid(True, linestyle = ':')

# Save the figure
fig11.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part4/4.2.svg', format='svg')
plt.savefig('A1 Aeroelastic design/Figures_part4/4.2.png', format='png')


'''FIGURE 4.3'''

DTU10MW_flex_path = './hawc_files/our_design/data/dtu_10mw_hawc2s_flex.pwr'
DTU10MW_flex_data = np.loadtxt(DTU10MW_flex_path, skiprows=1)

class DTU10MW_f_opt:
    V_o = DTU10MW_flex_data[:, 0]
    P = DTU10MW_flex_data[:, 1]
    T = DTU10MW_flex_data[:, 2]


# Set up a 1x2 grid layout
fig12, axes12 = plt.subplots(1, 2, figsize=(12, 4), dpi=500)

axes12[0].plot(DTU10MW_f_opt.V_o, DTU10MW_f_opt.P, color = colors[0], label='DTU 10MW (flexible) rotor')
axes12[0].plot(f_opt.V_o, f_opt.P, color = colors[1], linestyle='--', label='New IIIB (flexible) design')
axes12[0].set_ylabel(r"Power [kW]")
axes12[0].set_xlabel(r"Wind speed [m/s]")
axes12[0].legend()
axes12[0].grid(True, linestyle = ':')

axes12[1].plot(DTU10MW_f_opt.V_o, DTU10MW_f_opt.T, color = colors[0], label='DTU 10MW (flexible) rotor')
axes12[1].plot(f_opt.V_o, f_opt.T, color = colors[1], linestyle='--', label='New IIIB (flexible) design')
axes12[1].set_ylabel(r"Thrust [kN]")
axes12[1].set_xlabel(r"Wind speed [m/s]")
axes12[1].legend()
axes12[1].grid(True, linestyle = ':')

# Save the figure
fig12.tight_layout()
plt.savefig('A1 Aeroelastic design/Figures_part4/4.3.svg', format='svg')
plt.savefig('A1 Aeroelastic design/Figures_part4/4.3.png', format='png')