# %% Import modules
import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from aero_design_functions import get_design_functions_1MW, single_point_design, get_design_functions

# We're going from a class 1a to 3b
V_R_1a = 11.4
R_1a = 89.15 - 2.8
I_1a = 0.18
I_3b = 0.16

#Step 1
# Define symbolic variables
R_new, V_new = sp.symbols('R_new V_new')

# Define the equations using sympy
eq1 = (((V_R_1a * (1 + 2 * I_1a)) / (V_new * (1 + 2 * I_3b)))**(2/3)) * R_1a - R_new
eq2 = ((R_1a / R_new)**(2/3)) * V_R_1a - V_new

# Solve the system of equations
solution = sp.solve([eq1, eq2], (R_new, V_new))[0]
R_3b = float(solution[0])
V_R_3b = float(solution[1])

# Display the solution
print(f"Solved values: {R_3b} and {V_R_3b}")

#Step 4: Upscale the absolute blade-thickness and chord length

#Assuming a linear upscaling based on radius:
SF = R_3b/R_1a # Scaling factor (SF)

# The absolute blade-thickness and chord-length for the DTU 10MW reference turbine can be obtained for the _ae.dat file

ae_path = "C:\LAC_git\LAC_Assignments_Git\hawc_files\dtu_10mw\data\DTU_10MW_RWT_ae.dat"
data = np.loadtxt(ae_path, skiprows=2, usecols=[0,1,2])

# Now that we have the data, extract the radius, blade thickness, and chord length
R_1a_list = data[:, 0]  # Radius
c_1a_list = data[:, 1]  # Chord length
t_1a_list = data[:, 2]*data[:, 1]/100  # Absolute blade thickness

# Upscale based on the scaling factor
R_3b_list = R_1a_list * SF
t_3b_list = t_1a_list * SF
c_3b_list = c_1a_list * SF

# Visualize the scaling
plt.figure(0, figsize=(6, 4), dpi = 500)
plt.plot(R_1a_list, t_1a_list, label='Original Blade Thickness')
plt.plot(R_3b_list, t_3b_list, label='Scaled Blade Thickness', linestyle='--')
#plt.plot(R_1a_list, c_1a_list, label='Original Chord Length')
#plt.plot(R_3b_list, c_3b_list, label='Scaled Chord Length', linestyle='--')
plt.xlabel('Radius (m)')
plt.ylabel('Blade Thickness (m)')
plt.legend()
plt.grid()
plt.savefig("Figure_3_upscaling.png")

#Step 5: We use 3 blades

#Step 6: Use the example with the new inputs

pc_path = 'C:\LAC_git\LAC_Assignments_Git\hawc_files\dtu_10mw\data\DTU_10MW_RWT_pc.dat'
tc_data = np.loadtxt(pc_path, skiprows=2, usecols=[0,1,2])
tc_241_data = tc_data[52-2:69-2, :]
tc_301_data = tc_data[158-2:175-2, :]
tc_360_data = tc_data[264-2:276-2, :]
tc_480_data = tc_data[370-2:379-2, :]

# List of tc datasets for convenience
tc_datasets = [tc_241_data, tc_301_data, tc_360_data, tc_480_data]
titles = ['tc_241_data', 'tc_301_data', 'tc_360_data', 'tc_480_data']

# Plot for each tc dataset
for i, tc_data in enumerate(tc_datasets):
    alpha = tc_data[:, 0]  # Alpha (angle of attack)
    c_l = tc_data[:, 1]    # Lift coefficient (C_l)
    c_d = tc_data[:, 2]    # Drag coefficient (C_d)
    c_l_c_d = c_l / c_d    # Lift-to-drag ratio (C_l/C_d)

    # Create a figure for each tc_data set
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(f"Plots for {titles[i]}")
    
    # Plot 1: C_l vs C_d
    axs[0].plot(c_d, c_l, 'o-', color='blue')
    axs[0].set_title(f'C_l vs C_d')
    axs[0].set_xlabel('C_d')
    axs[0].set_ylabel('C_l')
    axs[0].grid()

    # Plot 2: C_l vs Alpha
    axs[1].plot(alpha, c_l, 'o-', color='green')
    axs[1].set_title(f'C_l vs Alpha')
    axs[1].set_xlabel('Alpha [deg]')
    axs[1].set_ylabel('C_l')
    axs[1].grid()

    # Plot 3: C_l vs C_l/C_d
    axs[2].plot(c_l_c_d, c_l, 'o-', color='red')
    axs[2].set_title(f'C_l vs C_l/C_d')
    axs[2].set_xlabel('C_l/C_d')
    axs[2].set_ylabel('C_l')
    axs[2].grid()

    # Display the plots
    plt.tight_layout()
    plt.show()

# %% Inputs
tsr = 9.0  # Tip-Speed-Ratio [-]
# A decrease in TSR will cause lower CT and CP and slight changes in planform
r_3b_hub = 2.8*SF  # Hub radius [m]
r_3b_list = np.linspace(r_3b_hub, R_3b - 0.1, 40)  # Rotor span [m]
c_3b_max = 6.2*SF  # Maximum chord size [m]
c_3b_root = 5.38*SF  # Chord size at the root [m]
B = 3  # Number of blades [#]
cl_scale = 1.0  # Change this value to scale the cl-values

# %% Plotting design functions
tc_plot = np.linspace(0, 100, 100)
cols = ['blue', 'orange', 'green']
fig1, axs1 = plt.subplots(3, 1, num=1)

cl_manu_points = [1.3, 1.23, 1.37, 0.54]
cd_manu_points = [0.013, 0.014, 0.021, 0.032]
alpha_manu_points = [8.1, 7.5, 5.8, 1.8]
tc_manu_points = [24.1, 30.1, 36, 48]

for i in range(1,4):
    cl_des_3b, cd_des_3b, aoa_des_3b, tc_vals_3b, cl_vals_3b, cd_vals_3b, aoa_vals_3b = get_design_functions(
        i
    )
    
    # %% Solving for the a single design
    chord, tc, twist, cl, cd, aoa, a, CLT, CLP, CT, CP = single_point_design(
        r_3b_list, t_3b_list, tsr, R_3b, cl_des_3b, cd_des_3b, aoa_des_3b, c_3b_root, c_3b_max, B
    )

    axs1[0].plot(tc_plot, cl_des_3b(tc_plot), color=cols[i-1], label = f"Design function {i}")
    axs1[0].scatter(tc_vals_3b, cl_vals_3b, color=cols[i-1],  marker='o')
    if i == 3:
        axs1[0].scatter(tc_manu_points, cl_manu_points, color='black',  marker='o')
    axs1[0].set_ylabel("$C_l$ [-]")
    axs1[0].set_xlim(0, 100)


    axs1[1].plot(tc_plot, cd_des_3b(tc_plot), color=cols[i-1])
    axs1[1].scatter(tc_vals_3b, cd_vals_3b, color=cols[i-1], marker='o')
    if i == 3:
        axs1[1].scatter(tc_manu_points, cd_manu_points, color='black',  marker='o')
    axs1[1].set_ylabel("$C_d$ [-]")
    axs1[1].set_xlim(0, 100)

    axs1[2].plot(tc_plot, aoa_des_3b(tc_plot), color=cols[i-1])
    axs1[2].scatter(tc_vals_3b, aoa_vals_3b, color=cols[i-1], marker='o')
    if i == 3:
        axs1[2].scatter(tc_manu_points, alpha_manu_points, color='black',  marker='o')
    axs1[2].set_ylabel(r"$\alpha$ [-]")
    axs1[2].set_xlabel(r"$t/c$ [deg]")
    axs1[2].set_xlim(0, 100)

fig1.tight_layout()
fig1.legend()
plt.show()