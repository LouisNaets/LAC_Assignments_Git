from lacbox.io import load_pwr, load_ind, load_inds
from lacbox.test import test_data_path
import matplotlib.pyplot as plt
import numpy as np

#Binary variables for including different parts of the script:
TSR_calc = 0

Step1_plot = 1      #plot 1wsp ind file: cl, l/d, AoA
Step2_1_plot = 0
Step2_2_plot = 0

if TSR_calc:
    TSRs = [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
    R = 88.93 #m
    v_o = 8 #m/s
    for TSR in TSRs:
        omega = TSR*v_o/R * (60/(2*np.pi))
        print(omega)

if Step1_plot:
    # Path for the file
    ind_path = "./hawc_files/group_7_design/res_hawc2s/group7_10mw_hawc2s_1wsp_u8000.ind"
    # Load the data
    ind_data = load_ind(ind_path)
    # Print the names in the dict
    print(ind_data.keys())
    #load design data
    design_data = np.load('3IIIB_design_data.npz',allow_pickle=True)
    # Access each array by its name
    r_des = design_data['r_des']
    cl_des = design_data['cl']
    cd_des = design_data['cd']
    aoa_des = design_data['aoa']

    fig, axs = plt.subplots(3, 1)

    # Cl
    axs[0].plot(ind_data["s_m"], ind_data["Cl"], label='HAWC2S $C_l$')
    axs[0].plot(r_des, cl_des, label='Design $C_l$')
    axs[0].set_ylabel("$C_l$ [-]")
    axs[0].legend()
    # CL/Cd
    axs[1].plot(ind_data["s_m"], ind_data["Cl"]/ind_data["Cd"])
    axs[1].plot(r_des, cl_des/cd_des)
    axs[1].set_ylabel("$C_l/C_d$ [-]")
    # aoa
    axs[2].plot(ind_data["s_m"], ind_data["aoa_rad"])
    axs[2].plot(r_des, aoa_des)
    axs[2].set_ylabel("Angle of attack [rad]")
    # Cl
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