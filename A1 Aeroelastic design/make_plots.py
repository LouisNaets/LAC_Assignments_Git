from lacbox.io import load_pwr, load_ind, load_inds
from lacbox.test import test_data_path
import matplotlib.pyplot as plt
import numpy as np

#Binary variables for including different parts of the script:
TSR_calc = 0

wsp1_plot = 0
Step2_1_plot = 1
Step2_2_plot = 0

if TSR_calc:
    TSR_range = np.arange(6, 10.1, 0.5)
    R = 88.93 #m
    v_o = 8 #m/s
    for TSR in TSRs:
        omega = TSR*v_o/R * (60/(2*np.pi))
        print(omega)


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

if wsp1_plot:
    # Path for the file
    ind_path = "./res_hawc2s/dtu_10mw_hawc2s_1wsp_u8000.ind"
    # Load the data
    ind_data = load_ind(ind_path)
    # Print the names in the dict
    print(ind_data.keys())

    fig, axs = plt.subplots(4, 2)

    # inflow angle
    axs[0, 0].plot(ind_data["s_m"], ind_data["flow_angle_rad"])
    axs[0, 0].set_ylabel("Inflow angle [rad]")
    # aoa
    axs[0, 1].plot(ind_data["s_m"], ind_data["aoa_rad"])
    axs[0, 1].set_ylabel("Angle of attack [rad]")
    # axial induction factor
    axs[1, 0].plot(ind_data["s_m"], ind_data["a"])
    axs[1, 0].set_ylabel("ax.-ind. ($a$) [-]")
    # tangential induction factor
    axs[1, 1].plot(ind_data["s_m"], ind_data["ap"])
    axs[1, 1].set_ylabel("tan.-ind. ($a_p$) [-]")
    # Cl
    axs[2, 0].plot(ind_data["s_m"], ind_data["Cl"])
    axs[2, 0].set_ylabel("$C_l$ [-]")
    # Cd
    axs[2, 1].plot(ind_data["s_m"], ind_data["Cd"])
    axs[2, 1].set_ylabel("$C_d$ [-]")
    # CP
    axs[3, 0].plot(ind_data["s_m"], ind_data["CP"])
    axs[3, 0].set_ylabel("Blade-span ($s$) [m]")
    axs[3, 0].set_ylabel("local-$C_P$ [-]")
    # CT
    axs[3, 1].plot(ind_data["s_m"], ind_data["CT"])
    axs[3, 1].set_ylabel("Blade-span ($s$) [m]")
    axs[3, 1].set_ylabel("local-$C_T$ [-]")

    fig.tight_layout()
    plt.show()

if Step2_1_plot:
    TSRs = [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
    # Path for the file
    pwr_path = "./res_hawc2s/dtu_10mw_hawc2s_multitsr.pwr"
    # Load the data
    pwr_data = load_pwr(pwr_path)

    # Path for the file
    inds_path = "./res_hawc2s/dtu_10mw_hawc2s_multitsr_u800%d.ind"
    # Load the data
    inds_data = load_inds([inds_path%i for i in range(7)])

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rc('legend', fontsize=14)
    plt.rc('axes', labelsize=16)
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)

    fig, axs = plt.subplots(4, 2, figsize=(6.4, 7), dpi=100)
    

    labels = ["%1.2f"%val for val in TSRs]

    # inflow angle
    axs[0, 0].plot(inds_data["s_m"], np.rad2deg(inds_data["flow_angle_rad"]), label=labels)
    axs[0, 0].set_ylabel(r"Inflow ($\phi$) [°]")
    # aoa
    axs[0, 1].plot(inds_data["s_m"], np.rad2deg(inds_data["aoa_rad"]), label=labels)
    axs[0, 1].set_ylabel(r"AoA ($\alpha$) [°]")
    # a
    axs[1, 0].plot(inds_data["s_m"], inds_data["a"], label=labels)
    axs[1, 0].set_ylabel("ax. ind. ($a$) [-]")
    # ap
    axs[1, 1].plot(inds_data["s_m"], inds_data["ap"], label=labels)
    axs[1, 1].set_ylabel("tan. ind. ($a_p$) [-]")
    # Cl
    axs[2, 0].plot(inds_data["s_m"], inds_data["Cl"], label=labels)
    axs[2, 0].set_ylabel("$C_l$ [-]")
    # Cd
    axs[2, 1].plot(inds_data["s_m"], inds_data["Cd"], label=labels)
    axs[2, 1].set_ylabel("$C_d$ [-]")
    axs[2, 1].legend(title=r"TSR ($\lambda$) [-]", ncol=2)
    # CP
    axs[3, 0].plot(inds_data["s_m"], inds_data["CP"], label=labels)
    axs[3, 0].set_ylabel("Blade-span ($s$) [m]")
    axs[3, 0].set_ylabel("local-$C_P$ [-]")
    axs[3, 0].set_xlabel("Blade-span ($s$) [m]")
    # CT
    axs[3, 1].plot(inds_data["s_m"], inds_data["CT"], label=labels)
    axs[3, 1].set_ylabel("Blade-span ($s$) [m]")
    axs[3, 1].set_ylabel("local-$C_T$ [-]")
    axs[3, 1].set_xlabel("Blade-span ($s$) [m]")

    fig.tight_layout()
    plt.show()

if Step2_2_plot:
    TSR = [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
    # Path for the file
    pwr_path = "./res_hawc2s/dtu_10mw_hawc2s_multitsr.pwr"
    # Load the data
    pwr_data = load_pwr(pwr_path)

    # Path for the file
    inds_path = "./res_hawc2s/dtu_10mw_hawc2s_multitsr_u800%d.ind"
    # Load the data
    inds_data = load_inds([inds_path%i for i in range(7)])

    plt.style.use("seaborn-v0_8-whitegrid")
    #print(plt.style.available)
    #plt.rcParams.update(plt.rcParamsDefault)

    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.4), dpi=200)


    # CP
    axs[0].plot(TSR, pwr_data["Cp"])
    axs[0].set_ylabel("$C_P$ [-]")
    axs[0].set_xlabel(r"TSR ($\lambda$) [-]")
    # CT
    axs[1].plot(TSR, pwr_data["Ct"])
    axs[1].set_ylabel("$C_T$ [-]")
    axs[1].set_xlabel(r"TSR ($\lambda$) [-]")

    fig.tight_layout()
    plt.show()