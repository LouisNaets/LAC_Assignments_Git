# -*- coding: utf-8 -*-
from lacbox.io import load_stats, load_oper
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import statistics as stats
from scipy.stats import weibull_min

# Set constants and turbine parameters
WTG = 'DTU10MW'  # Turbine model
N_T = 630720000  # Total cycles over turbine lifetime
n_eq = 10e6  # Equivalent cycles for 10-minute DEL
U_ave = 10  # Average wind speed for Weibull distribution
U_std = 2   # Wind speed standard deviation
c = 2 / np.sqrt(np.pi) * U_ave  # Scale parameter for Weibull
k = 2  # Shape parameter for Weibull

# For verification
extreme_design_loads = {
    'TbFA': 364297.47 / 1.35,
    'TbSS': 117720.32 / 1.35,
    'YbTilt': 54987.89 / 1.35,
    'YbRoll': 23907.46 / 1.35,
    'ShftTrs': -20365.17 * 0.7,
    'OoPBRM': -72275.14,
    'IPBRM': 40770.71
}

fatigue_design_loads = {
    'TbFA': 129826.23,
    'TbSS': 55207.15,
    'YbTilt': 32488.02,
    'YbRoll': 4039.54,
    'ShftTrs': 2839.47,
    'OoPBRM': 31335.80,
    'IPBRM': 31017.46
}


m_values = {
    'TbFA': 4, 'TbSS': 4, 'YbTilt': 4, 'YbRoll': 4,
    'ShftTrs': 4, 'OoPBRM': 10, 'IPBRM': 10
}

if WTG == 'group7':
    # analysis settings
    HAWC2S_PATH = './hawc_files/our_design/data/group7_3B_design_flex.opt'  # path to .pwr or .opt file
    STATS_PATH = './A4 Design Loads and AEP/Assignment/group7_turbB_stats.csv'  # path to mean steady stats
    SUBFOLDER = 'tcb'  # which subfolder to plot: tca or tcb

    # what channels we want to get the design loads
    chan_ids = ['TbFA', 'TbSS', 'YbTilt', 'YbRoll', 'ShftTrs', 'OoPBRM', 'IPBRM']

    turb_ids = ['path', 'filename', 'subfolder', 'ichan', 'names', 'units', 'desc',
                'mean', 'max', 'min', 'std', '1%', '50%', '99%', 'del3', 'del4', 'del5',
                'del8', 'del10', 'del12', 'wsp']

    # load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
    df, wsps = load_stats(STATS_PATH, statstype='turb')

if WTG == 'DTU10MW':
    # analysis settings
    HAWC2S_PATH = './hawc_files\dtu_10mw\data\dtu_10mw_flex.opt'  # path to .pwr or .opt file
    STATS_PATH = './A4 Design Loads and AEP\Assignment\dtu_10mw_turb_stats.hdf5'  # path to mean steady stats
    SUBFOLDER = 'tcb'  # which subfolder to plot: tca or tcb

    # what channels we want to get the design loads
    chan_ids = ['TbFA', 'TbSS', 'YbTilt', 'YbRoll', 'ShftTrs', 'OoPBRM', 'IPBRM']

    turb_ids = ['path', 'filename', 'subfolder', 'ichan', 'names', 'units', 'desc',
                'mean', 'max', 'min', 'std', '1%', '50%', '99%', 'del3', 'del4', 'del5',
                'del8', 'del10', 'del12', 'wsp']

    # load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
    df, wsps = load_stats(STATS_PATH, subfolder = SUBFOLDER, statstype='turb')

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

# load/calc. the stuff we need from the HAWC2S opt/pwr file for the operational data comparisons
opt_dict = load_oper(HAWC2S_PATH)
h2s_u, h2s_pitch, h2s_rotspd, = opt_dict['ws_ms'], opt_dict['pitch_deg'], opt_dict['rotor_speed_rpm']
h2s_paero, h2s_thrust = opt_dict['power_kw'], opt_dict['thrust_kn']
h2s_aerotrq = h2s_paero / (h2s_rotspd * np.pi / 30)

# get hawc2 thrust and aerodynamic torque for theoretical calculations
h2_thrust = df.filter_channel('Thrust', CHAN_DESCS)['mean']
h2_aero_trq = df.filter_channel('GenTrq', CHAN_DESCS)['mean'] / GENEFF * 1e-3  # aerodynamic torque [kNm]

# Weibull bin probabilities
bin_edges = np.array([4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 
                      14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5])
bin_probabilities = np.array([
    weibull_min.cdf(bin_edges[i + 1], k, scale=c) - weibull_min.cdf(bin_edges[i], k, scale=c)
    for i in range(len(bin_edges) - 1)
])
bin_probabilities /= np.sum(bin_probabilities)
n_i = N_T * bin_probabilities  # Lifetime cycles for each bin

# Helper function to calculate 10-minute combined DEL for each bin
def combine_10min_DELs(DELs, m):
    return (np.mean(DELs ** m)) ** (1 / m)

# Helper function to calculate the lifetime fatigue load
def calculate_lifetime_fatigue(DELs_per_bin, n_i, m):
    adjusted_DELs = (n_i / n_eq) * (DELs_per_bin ** m)
    sum_adjusted_DELs = np.sum(adjusted_DELs)
    return sum_adjusted_DELs ** (1 / m)

# Loop over each channel and plot the steady state with the theory line
for iplot, chan_id in enumerate(chan_ids):
    
    # Isolate the channel data
    chan_df = df.filter_channel(chan_id, CHAN_DESCS)
    Wohler_exponent = m_values[chan_id]

    # Extract HAWC2 wind and the stats ('mean', 'min', 'max') for the channel
    h2_wind = np.array(chan_df['wsp'])
    HAWC2val_mean = np.array(chan_df['mean'])
    HAWC2val_min = np.array(chan_df['min'])
    HAWC2val_max = np.array(chan_df['max'])
    HAWC2val_del = np.array(chan_df[f'del{m_values[chan_id]}'])

    # Sort HAWC2 values by increasing wind speed
    i_h2 = np.argsort(h2_wind)
    h2_wind_sorted = h2_wind[i_h2]
    HAWC2val_mean_sorted = HAWC2val_mean[i_h2]
    HAWC2val_min_sorted = HAWC2val_min[i_h2]
    HAWC2val_max_sorted = HAWC2val_max[i_h2]
    HAWC2val_del_sorted = np.array([HAWC2val_del[i:i+6] for i in range(0, len(HAWC2val_del), 6)])

    # Calculate the max of the max, and min of min
    unique_wind_speeds = np.unique(h2_wind_sorted)
    mean_HAWC2val_per_wind = [np.mean(HAWC2val_mean_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
    min_HAWC2val_per_wind = [np.min(HAWC2val_min_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
    max_HAWC2val_per_wind = [np.max(HAWC2val_max_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]

    # Calculate characteristic values including safety factors
    SF = 1.35  # safety factor
    PSF = 1.25  # partial safety factor
    if chan_id == 'OoPBRM' or chan_id == 'IPBRM':
        max_charval_per_wind = SF * PSF * np.array(max_HAWC2val_per_wind)
        min_charval_per_wind = SF * PSF * np.array(min_HAWC2val_per_wind)
    
    else:
        max_charval_per_wind = PSF * np.array(max_HAWC2val_per_wind)
        min_charval_per_wind = PSF * np.array(min_HAWC2val_per_wind)

    # Store the max and min values for the characteristic values for the channel
    max_charval = np.max(max_charval_per_wind)
    min_charval = np.min(min_charval_per_wind)
    abs_max_charval = max(abs(max_charval), abs(min_charval))

    print(f'{chan_id}: Ultimate Design Load = {abs_max_charval:.2f} kNm')

    # Plot the max values max_HAWC2val_per_wind vs wind speed (unique_wind_speeds) as bars
    plt.figure()
    plt.bar(unique_wind_speeds, max_charval_per_wind, width=0.5, color='b', alpha=0.5, label='HAWC2 max')
    plt.bar(unique_wind_speeds, min_charval_per_wind, width=0.5, color='r', alpha=0.5, label='HAWC2 min')
    plt.axhline(extreme_design_loads[chan_id], color='r', linestyle='--', label='DTU 10MW Extreme Design Load')
    plt.title(chan_id)
    plt.xlabel('Wind speed [m/s]')
    plt.ylabel('Max value [kNm]')
    plt.grid()
    plt.legend()

    # Calculate combined 10-min DEL and lifetime fatigue load
    DELs_10min = np.array([combine_10min_DELs(DEL_bin, Wohler_exponent) for DEL_bin in HAWC2val_del_sorted])
    lifetime_fatigue_load = calculate_lifetime_fatigue(DELs_10min, n_i, Wohler_exponent)

    # print(f"{chan_id}: 10-min Combined DELs:")
    # for i, DEL in enumerate(DELs_10min):
        # print(f"  Wind speed bin {h2_wind_sorted[i*6]}: {DEL:.8f} kNm")

    print(f"{chan_id}: Lifetime Fatigue Load = {lifetime_fatigue_load:.8f} kNm")


plt.show()



    