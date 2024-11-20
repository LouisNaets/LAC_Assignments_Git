# -*- coding: utf-8 -*-
"""Compare mean value of steady simulations in HAWC2 (blue dots) to "theory" (lines).

For the operational parameters (i.e., pitch, rotor speed, etc.), the "theoretical" values
are the corresponding values in a .pwr file.

For the load channels of interest, the theoretical lines are calculated according to the
theoretical equations we derived for each load channel, as a function of thrust/torque/
gravity moment.

YOUR TASK! Add the lines that calculate the theory, as prompted by the slides.
"""
from lacbox.io import load_stats, load_oper
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import statistics as stats

def AEP_data(STATS_PATH, HAWC2S_PATH, SUBFOLDER, ubar, winddist):
    '''
    Calculate the AEP of the turbine based on the power output and the wind speed distribution
    '''
    # load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
    # turbine constants
    GENEFF = 0.94  # generator/gearbox efficienty [%]
    FG_TIMES_DY = 6250  # yaw-bearing pitch moment due to gravity [kNm]
    if 'notilt' in SUBFOLDER:
        DZ_YB = 2.75  # distance from hub center to yaw bearing [m]
        DZ_TB = 115.63 + DZ_YB  # distance from hub center to tower base [m]
    else:
        DZ_YB = 2.75 + 7.1*np.sin(5*np.pi/180)  # distance from hub center to yaw bearing [m]
        DZ_TB = 115.63 + DZ_YB  # distance from hub center to tower base [m]
    CHAN_DESCS = {'BldPit': 'pitch1 angle',  # dictionary used to identify which descriptions
                'RotSpd': 'rotor speed',  # in the HAWC2 statistics file correspond to which
                'Thrust': 'aero rotor thrust',  # channels we want
                'GenTrq': 'generator torque',
                'ElPow': 'pelec',
                'TbFA': 'momentmx mbdy:tower nodenr:   1',
                'TbSS': 'momentmy mbdy:tower nodenr:   1',
                'YbTilt': 'momentmx mbdy:tower nodenr:  11',
                'YbRoll': 'momentmy mbdy:tower nodenr:  11',
                'ShftTrs': 'momentmz mbdy:shaft nodenr:   4',
                'OoPBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: hub1',
                'IPBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: hub1',
                'FlpBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: blade1',
                'EdgBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: blade1',
                'OoPHub': 'momentmx mbdy:hub1 nodenr:   1 coo: hub1',
                'IPHub': 'momentmy mbdy:hub1 nodenr:   1 coo: hub1',
                }

    # what channels we want to plot
    chan_ids = ['BldPit', 'RotSpd', 'Thrust', 'GenTrq', 'ElPow', 'TbFA', 'TbSS',
                'YbTilt', 'YbRoll', 'ShftTrs', 'OoPBRM', 'IPBRM']

    turb_ids = ['path', 'filename', 'subfolder', 'ichan', 'names', 'units', 'desc',
                'mean', 'max', 'min', 'std', '1%', '50%', '99%', 'del3', 'del4', 'del5',
                'del8', 'del10', 'del12', 'wsp']
    if SUBFOLDER=='':
        df, wsps = load_stats(STATS_PATH, statstype='turb')
    else:
        df, wsps = load_stats(STATS_PATH, subfolder=SUBFOLDER, statstype='turb')

    # load/calc. the stuff we need from the HAWC2S opt/pwr file for the operational data comparisons
    opt_dict = load_oper(HAWC2S_PATH)
    h2s_u, h2s_pitch, h2s_rotspd, = opt_dict['ws_ms'], opt_dict['pitch_deg'], opt_dict['rotor_speed_rpm']
    h2s_paero, h2s_thrust = opt_dict['power_kw'], opt_dict['thrust_kn']
    h2s_aerotrq = h2s_paero / (h2s_rotspd * np.pi / 30)

    # get hawc2 thrust and aerodynamic torque for theoretical calculations
    h2_thrust = df.filter_channel('Thrust', CHAN_DESCS)['mean']
    h2_aero_trq = df.filter_channel('GenTrq', CHAN_DESCS)['mean'] / GENEFF * 1e-3  # aerodynamic torque [kNm]

    # ---------------- AEP calculations --------------------
    v = np.arange(5,25)
    bin_edges = np.arange(4.5, 25.5)

    if winddist == 'weibull':
        TI = 0.18 #TI = sigma/mean
        #sigma = TI * ubar #standard deviation
        k = 2                                   # # Weibull shape parameter
        c = 2/np.sqrt(np.pi)*ubar             # Weibull scale parameter
        bin_prop =  (k / c) * (bin_edges/ c)**(k - 1) * np.exp(-(bin_edges / c)**k)
        bin_prop = np.exp(-((v-0.5)/(2*ubar/np.sqrt(np.pi)))**2) - np.exp(-((v+0.5)/(2*ubar/np.sqrt(np.pi)))**2) #weibull equation SLIDE 14
        print('WEIBULL distribution for bin probaility selected')
        print("Weibull scale parameter (c):", round(c, 3))
        print("Weibull shape parameter (k):", round(k, 3))
    
    elif winddist == 'rayleigh':
        bin_prop = 1-np.exp(-np.pi*(v/(2*ubar))**2)          #IEC rayleigh equation 8
        bin_prop = np.pi/2*(v/ubar**2)*np.exp((-np.pi/4)*(v/ubar)**2)      #paper  Rayleigh pdf
        print('RAYLEIGH distribution for bin probaility selected')
    else:   
        print('ERROR: Invalid wind distribution selected')
        exit()

    bin_prop_jenni =  [0.06442809, 0.0709227,  0.07472189, 0.07591763, 0.0747455,  0.07155162,
    0.06675382, 0.06080169, 0.05413969, 0.04717648, 0.04026251, 0.03367663,
    0.02762128, 0.02222493, 0.01755019, 0.01360521, 0.01035688, 0.00774379,
    0.00568809, 0.0041053 ]

    # bin_prop_avg = np.array([])
    # for i in range(0, len(bin_prop)-1):
    #     np.append(bin_prop_avg, np.mean(bin_prop[i:i+1]))
    # bin_prop = bin_prop_avg
    print("Wind speeds: ", v)
    print("Bin probabilties: ",bin_prop)
    print('Sum of bin probability: ', sum(bin_prop))
    #print('Weibull jenni sum: ', sum(bin_prop_jenni))
    print("Element-wise difference: ", bin_prop-np.array(bin_prop_jenni))

    # plt.plot(v, weibull, label='Weibull distribution')
    # plt.plot(v, bin_prop_jenni, label='Bin prop Jenni')
    # plt.legend()
    #plt.show()

    #get hawc2 electric power
    h2_power = df.filter_channel('ElPow', CHAN_DESCS)
    h2_power_mean = np.array(h2_power['mean'])
    # Extract HAWC2 wind and the stats ('mean', 'min', 'max') for the channel
    h2_wind = np.array(h2_power['wsp'])
    i_h2 = np.argsort(h2_wind)
    h2_wind_sorted = h2_wind[i_h2]
    h2_power_mean_sorted = h2_power_mean[i_h2]

    #Average electrical power over turbulent seeds
    h2_power_mean_turb = np.array([])
    for i in range(0, len(h2_power_mean_sorted), 6):  
        mean = stats.mean(h2_power_mean_sorted[i:i+6])
        h2_power_mean_turb = np.append(h2_power_mean_turb, mean)  
    print('Electric power:', h2_power_mean_turb/1e6)
    p_tot = h2_power_mean_turb*bin_prop
    R = 1
    P_tot_reliability = np.sum(p_tot) * R
    AEP = P_tot_reliability * 8760*1e-9
    print('AEP:', str(AEP), 'GWh')
    print('')
    print('----------------------End of AEP calculation----------------------')
    print('')
    # plt.figure(figsize=(6, 4), clear=True)
    # plt.plot(v, p_tot*1e-9, label = 'Mean')
    # plt.legend(fontsize=18)
    # plt.xlabel('Wind speed [m/s]', fontsize=18)
    # plt.ylabel('Power [MWh]' ,fontsize=18)
    # plt.xticks(fontsize=18)
    # plt.yticks(fontsize=18)
    # plt.show()

    return AEP, p_tot, bin_prop, h2_wind_sorted, v, h2_power_mean_turb

# analysis settings
HAWC2S_PATH_DTU = './A4 Design Loads and AEP/Assignment/dtu_10mw_flex_minrotspd.opt'  # path to .pwr or .opt file
STATS_PATH_DTU = './A4 Design Loads and AEP/Assignment/dtu_10mw_turb_stats.hdf5'  # path to mean steady stats
SUBFOLDER_DTU = 'tcb'  # which subfolder to plot: tca or tcb

HAWC2S_PATH_OURS = './hawc_files/our_design/data/group7_3B_design_flex.opt'  # path to .pwr or .opt file
STATS_PATH_OURS = './A4 Design Loads and AEP/Assignment/group7_turbB_stats_ts.csv'  # path to mean steady stats
SUBFOLDER_OURS = ''  # which subfolder to plot: tca or tcb


AEP_DTU, p_tot_DTU, bin_prop_DTU, h2_wind_sorted_DTU, v_bins_DTU, power_mean_DTU = AEP_data(STATS_PATH_DTU, HAWC2S_PATH_DTU, SUBFOLDER_DTU, 7.5, 'weibull')
AEP_OURS, p_tot_OURS, bin_prop_OURS, h2_wind_sorted_OURS, v_bins_OURS, power_mean_OURS = AEP_data(STATS_PATH_OURS, HAWC2S_PATH_OURS, SUBFOLDER_OURS, 7.5, 'weibull')


# --- Power Production per Wind Bin (Bar Plot) ---

# Interpolated power curves
total_power_curve_DTU = np.interp(v_bins_DTU, v_bins_DTU, p_tot_DTU)
total_power_curve_OURS = np.interp(v_bins_OURS, v_bins_OURS, p_tot_OURS)

# Initialize the figure
plt.figure(figsize=(12, 8))

# --- Power Production per Wind Bin (Bar Plot) ---
# DTU 10 MW - Dark Blue
plt.bar(v_bins_DTU, p_tot_DTU*1e-6, width=0.8, alpha=0.7, color='navy', label='Power Output - DTU 10 MW')
# Custom Design - Light Blue
plt.bar(v_bins_OURS, p_tot_OURS*1e-6, width=0.8, alpha=0.7, color='skyblue', label='Power Output - Custom Design')

# --- Total Power Curve (Line Plot) ---
# DTU 10 MW - Solid Line
plt.plot(v_bins_DTU, power_mean_DTU*1e-6, 'o-', color='darkblue', label='Total Power Curve - DTU 10 MW')
# Custom Design - Dashed Line
plt.plot(v_bins_OURS, power_mean_OURS*1e-6, 'o--', color='deepskyblue', label='Total Power Curve - Custom Design')

# --- Bin Probabilities (Overlaid Bar Plot) ---
# DTU 10 MW - Dark Blue (transparent)
plt.bar(v_bins_DTU, bin_prop_DTU*100, width=0.4, alpha=0.3, color='navy', label='Bin Probabilities - DTU 10 MW')
# Custom Design - Light Blue (transparent)
plt.bar(v_bins_OURS, bin_prop_OURS*100, width=0.4, alpha=0.3, color='skyblue', label='Bin Probabilities - Custom Design')

# --- Labels, Legends, and Formatting ---
plt.xlabel('Wind Speed [m/s]')
plt.ylabel('Power [MW] / Probability')
plt.title('Comparison of Power Production and Bin Probabilities')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Display the plot
#plt.show()



# Initialize the figure and create two y-axes
fig, ax1 = plt.subplots(figsize=(12, 8))

# --- Total Power Curve (Line Plot) on the left y-axis (ax1) ---
ax1.plot(v_bins_DTU, power_mean_DTU* 1e-6, 'o-', color='indigo', label='Total Power Curve - DTU 10 MW')
ax1.plot(v_bins_OURS, power_mean_OURS * 1e-6, 'o--', color='mediumpurple', label='Total Power Curve - Our Design')

# Label and formatting for the left y-axis
ax1.set_xlabel('Wind Speed [m/s]', fontsize=18)
ax1.set_ylabel('Power Output [MW] & Bin probability [%]', color='darkorchid', fontsize=18)
ax1.tick_params(axis='y', labelcolor='darkorchid', labelsize=18)
ax1.tick_params(axis='x', labelsize=18)

ax1.plot(v_bins_DTU, bin_prop_DTU*100, 's-', color='indigo', label='Bin Probabilities - DTU 10 MW')
ax1.plot(v_bins_OURS, bin_prop_OURS*100, 's--',color='mediumpurple', label='Bin Probabilities - Our design')
#ax1.set_ylabel('Bin Probability', color='purple', fontsize=18)

# Create a secondary y-axis for weighted power (ax2)
ax2 = ax1.twinx()

# --- Power Production per Wind Bin (Bar Plot) on the left y-axis (ax1) ---
ax2.bar(v_bins_DTU, p_tot_DTU * 8670 * 1e-6, width=0.8, alpha=0.7, color='navy', label='Power Output - DTU 10 MW')
ax2.bar(v_bins_OURS, p_tot_OURS * 8670 * 1e-6, width=0.8, alpha=0.7, color='skyblue', label='Power Output - Our Design')

# Label and formatting for the right y-axis
ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=18)
ax2.set_ylabel('Weighted Power Output [MWh]', color='darkblue', fontsize=18)

# Legends for both axes
ax1.legend(loc='center right', bbox_to_anchor=(1, 0.7), fontsize=16)
ax2.legend(loc='center right', bbox_to_anchor=(1, 0.55), fontsize=16)

# Title and grid
#plt.title('Comparison of Power Production and Bin Probabilities')
ax1.grid(True)
plt.tight_layout()

# Display the plot
plt.savefig('./A4 Design Loads and AEP/Assignment/Figures/AEP_comparison.png', dpi=300)
plt.savefig('./A4 Design Loads and AEP/Assignment/Figures/AEP_comparison.svg', dpi=300)
plt.show()










# '''AEP Q1 - Weibull parameters'''
# U_mean =10
# sigma_U = 5

# C = (2/np.sqrt(np.pi))*U_mean
# k = (sigma_U/U_mean)**(-1.086)
# print("Weibull scale parameter (C):", round(C, 3))
# print("Weibull shape parameter (k):", round(k, 3))
# print('Note: These expressions are valid only when 1.6 < k < 3')

# '''AEP Q2 - Simple power curve AEP'''
# V_bins = [0,5,8,12,14,25,100] #100 is just a high random number
# P_bins = [0,20,35,40,45,0]

# V_prob = []
# for i in range(0,len(V_bins)-1):
#     V_prob.append(np.exp(-(V_bins[i]/C)**k) - np.exp(-(V_bins[i+1]/C)**k))

# P = [p * v for p, v in zip(P_bins, V_prob)]
# P_tot = sum(P)

# print("Power generation before reliability:", round(P_tot,3), 'kW')

# R = 0.95
# P_tot_reliability = P_tot * R
# AEP = P_tot_reliability * 365.25 * 24/1000

# print('AEP:', str(round(AEP,2)), 'MWh')






# ---------------------- Plotting ----------------------
# initialize the figure and axes
# fig, axs = plt.subplots(4, 3, figsize=(12, 8), clear=True)

# # Set the opacity and marker size variables
# dot_opacity = 0.25  # Opacity for individual points
# line_opacity = 0.8  # Opacity for the mean lines
# dot_size = 10       # Size for individual points
# line_size = 40      # Size for mean points
# max_color = 'tab:gray'
# mean_color = 'tab:blue'
# min_color = 'tab:orange'

# # Loop over each channel and plot the steady state with the theory line
# for iplot, chan_id in enumerate(chan_ids):
    
#     # Isolate the channel data
#     chan_df = df.filter_channel(chan_id, CHAN_DESCS)

#     # Extract HAWC2 wind and the stats ('mean', 'min', 'max') for the channel
#     h2_wind = np.array(chan_df['wsp'])
#     HAWC2val_mean = np.array(chan_df['mean'])
#     HAWC2val_min = np.array(chan_df['min'])
#     HAWC2val_max = np.array(chan_df['max'])

#     # Sort HAWC2 values by increasing wind speed
#     i_h2 = np.argsort(h2_wind)
#     h2_wind_sorted = h2_wind[i_h2]
#     HAWC2val_mean_sorted = HAWC2val_mean[i_h2]
#     HAWC2val_min_sorted = HAWC2val_min[i_h2]
#     HAWC2val_max_sorted = HAWC2val_max[i_h2]

#     # Calculate the mean for each unique wind speed for 'mean', 'min', and 'max'
#     unique_wind_speeds = np.unique(h2_wind_sorted)
#     mean_HAWC2val_per_wind = [np.mean(HAWC2val_mean_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
#     min_HAWC2val_per_wind = [np.mean(HAWC2val_min_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
#     max_HAWC2val_per_wind = [np.mean(HAWC2val_max_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]

#     # Plot the results
#     ax = axs.flatten()[iplot]
#     # Individual 'mean' values with lower opacity and smaller markers
#     ax.scatter(h2_wind_sorted, HAWC2val_mean_sorted, color=mean_color, alpha=dot_opacity, s=dot_size, label='Mean')
#     # Mean of 'mean' with higher opacity and larger markers
#     ax.plot(unique_wind_speeds, mean_HAWC2val_per_wind, color=mean_color, alpha=line_opacity, markersize=line_size)#, label='Mean of means')

#     # Individual 'min' values
#     ax.scatter(h2_wind_sorted, HAWC2val_min_sorted, color=min_color, alpha=dot_opacity, s=dot_size, label='Min')
#     # Mean of 'min'
#     ax.plot(unique_wind_speeds, min_HAWC2val_per_wind, color=min_color, alpha=line_opacity, markersize=line_size)#, label='Mean of min')

#     # Individual 'max' values
#     ax.scatter(h2_wind_sorted, HAWC2val_max_sorted, color=max_color, alpha=dot_opacity, s=dot_size, label='Max')
#     # Mean of 'max'
#     ax.plot(unique_wind_speeds, max_HAWC2val_per_wind, color=max_color, alpha=line_opacity, markersize=line_size)#, label='Mean of max')

#     # Formatting the plot
#     ax.grid('on')
#     ax.set(xlabel='Wind speed [m/s]' if iplot > 8 else None,
#            ylabel=f'{chan_id} [{chan_df.units.iloc[0]}]', xlim=[4, 25])

# # Add legends and format the figure
# axs[0, 0].legend()
# #axs[1, 2].legend()
# fig.suptitle(f'Case: DTU 10MW turbine - {SUBFOLDER}')
# fig.tight_layout()
# plt.show()



