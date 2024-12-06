# -*- coding: utf-8 -*-
from lacbox.io import load_stats, load_oper
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import statistics as stats
from scipy.stats import weibull_min

# Set constants and turbine parameters
WTG = 'Redesign'  # Turbine model
N_T = 630720000  # Total cycles over turbine lifetime
n_eq = 10e6  # Equivalent cycles for 10-minute DEL

# For verification
extreme_design_loads = {
    'TbFA': 364297.47 / 1.35,
    'TbSS': 117720.32 / 1.35,
    'YbTilt': 54987.89 / 1.35,
    'YbRoll': 23907.46 / 1.35,
    'ShftTrs': -20365.17 / 1.35,
    'OoPBRM': -72275.14,
    'IPBRM': 40770.71,
    'FlpBRM': -53583.14,
    'EdgBRM': 21464.92
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
    'ShftTrs': 4, 'OoPBRM': 10, 'IPBRM': 10, 'EdgBRM': 10, 'FlpBRM': 10
}

if WTG == 'Redesign':
    U_ave = 7.5  # Average wind speed for Weibull distribution
    # analysis settings
    HAWC2S_PATH = './hawc_files/individual_design/data/individual_design_flex_minrotspd.opt'  # path to .pwr or .opt file
    STATS_PATH = './A5 Individual design/Part 4/individual_design_turb_tcb_stats.csv'  # path to mean steady stats
    SUBFOLDER = 'tcb'  # Turbulence Class B

    # load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
    df, wsps = load_stats(STATS_PATH, statstype='turb')

if WTG == 'DTU10MW':
    U_ave = 10  # Average wind speed for Weibull distribution
    # analysis settings
    HAWC2S_PATH = './hawc_files/dtu_10mw/data/dtu_10mw_flex_minrotspd.opt'  # path to .pwr or .opt file
    STATS_PATH = './A5 Individual design/Part 4/dtu_10mw_turb_stats.hdf5'  # path to mean steady stats
    SUBFOLDER = 'tca'  # Turbulence Class A

    # load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
    df, wsps = load_stats(STATS_PATH, subfolder = SUBFOLDER, statstype='turb')

# what channels we want to get the design loads
chan_ids = ['TbFA', 'TbSS', 'YbTilt', 'YbRoll', 'ShftTrs', 'OoPBRM', 'IPBRM', 'EdgBRM', 'FlpBRM']

turb_ids = ['path', 'filename', 'subfolder', 'ichan', 'names', 'units', 'desc',
            'mean', 'max', 'min', 'std', '1%', '50%', '99%', 'del3', 'del4', 'del5',
            'del8', 'del10', 'del12', 'wsp']

CHAN_DESCS = {'TbFA': 'momentmx mbdy:tower nodenr:   1',
              'TbSS': 'momentmy mbdy:tower nodenr:   1',
              'YbTilt': 'momentmx mbdy:tower nodenr:  11',
              'YbRoll': 'momentmy mbdy:tower nodenr:  11',
              'ShftTrs': 'momentmz mbdy:shaft nodenr:   4',
              'OoPBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: hub1',
              'IPBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: hub1',
              'FlpBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: blade1',
              'EdgBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: blade1',
              'OoPHub': 'momentmx mbdy:hub1 nodenr:   1 coo: hub1',
              'IPHub': 'momentmy mbdy:hub1 nodenr:   1 coo: hub1'
              #'EdgBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped',
              #'FlpBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped'
              }

# Weibull bin probabilities
c = 2 / np.sqrt(np.pi) * U_ave  # Scale parameter for Weibull
k = 2  # Shape parameter for Weibull

bin_edges = np.array([4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 
                      14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5])
bin_probabilities = np.array([
    weibull_min.cdf(bin_edges[i + 1], k, scale=c) - weibull_min.cdf(bin_edges[i], k, scale=c)
    for i in range(len(bin_edges) - 1)
])

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
fig, axs = plt.subplots(3, 3, figsize=(15, 15))
axs = axs.flatten()

for iplot, chan_id in enumerate(chan_ids):
    
    # Isolate the channel data
    chan_df = df.filter_channel(chan_id, CHAN_DESCS)
    Wohler_exponent = m_values[chan_id]

    # Extract HAWC2 wind and the stats ('mean', 'min', 'max') for the channel
    h2_wind = np.array(chan_df['wsp'])
    HAWC2val_mean = np.array(chan_df['mean'])
    HAWC2val_min = np.array(chan_df['min'])
    HAWC2val_max = np.array(chan_df['max'])
    HAWC2val_del = np.array(chan_df[f'del{Wohler_exponent}'])
    # print(f'{chan_id}: {(h2_wind)}')
    # print(f'{chan_id}: {(HAWC2val_max)}')

    # Sort HAWC2 values by increasing wind speed
    i_h2 = np.argsort(h2_wind)
    h2_wind_sorted = h2_wind[i_h2]
    HAWC2val_mean_sorted = HAWC2val_mean[i_h2]
    HAWC2val_min_sorted = HAWC2val_min[i_h2]
    HAWC2val_max_sorted = HAWC2val_max[i_h2]
    HAWC2val_del_sorted = HAWC2val_del[i_h2]
    
    # Calculate the mean of the max, and mean of min
    unique_wind_speeds = np.unique(h2_wind_sorted)
    mean_HAWC2val_per_wind = [np.mean(HAWC2val_mean_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
    min_HAWC2val_per_wind = [np.mean(HAWC2val_min_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
    max_HAWC2val_per_wind = [np.mean(HAWC2val_max_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
    HAWC2val_del_per_wind = np.array([HAWC2val_del_sorted[h2_wind_sorted == wsp] for wsp in unique_wind_speeds])
    # print(f'{chan_id}: {np.max(max_HAWC2val_per_wind)}')

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
    ax = axs[iplot]
    for wsp, max_val, min_val in zip(unique_wind_speeds, max_charval_per_wind, min_charval_per_wind):
        if max_val >= min_val:
            ax.bar(wsp, max_val, width=0.5, color='navy', alpha=0.8, label='HAWC2 max' if wsp == unique_wind_speeds[0] else "")
            ax.bar(wsp, min_val, width=0.5, color='skyblue', alpha=0.8, label='HAWC2 min' if wsp == unique_wind_speeds[0] else "")
        else:
            ax.bar(wsp, max_val, width=0.5, color='navy', alpha=0.8, label='HAWC2 max' if wsp == unique_wind_speeds[0] else "")
            ax.bar(wsp, min_val, width=0.5, color='skyblue', alpha=0.8, label='HAWC2 min' if wsp == unique_wind_speeds[0] else "")
    ax.axhline(extreme_design_loads[chan_id], color='tab:gray', linestyle='--', label='DTU 10MW Extreme Design Load')
    ax.set_title(chan_id)
    if iplot >= 6:
        ax.set_xlabel('Wind speed [m/s]')
    if iplot % 2 == 0:
        ax.set_ylabel('Max value [kNm]')
    ax.grid()
    if iplot == 0:
        ax.legend()

    # Calculate combined 10-min DEL and lifetime fatigue load
    # Calculate 10-minute DELs for each bin
    DELs_10min = []

    for DEL_bin in HAWC2val_del_per_wind:
        combined_DEL = combine_10min_DELs(DEL_bin, Wohler_exponent)
        DELs_10min.append(combined_DEL)

    DELs_10min = np.array(DELs_10min)
    lifetime_fatigue_load = calculate_lifetime_fatigue(DELs_10min, n_i, Wohler_exponent)

    print(f"{chan_id}: Lifetime Fatigue Load = {lifetime_fatigue_load:.8f} kNm")

fig.tight_layout()
plt.savefig(f'./A5 Individual design/Figures 4/{WTG} Ultimate design_loads.png')
plt.savefig(f'./A5 Individual design/Figures 4/{WTG} Ultimate design_loads.svg')
# plt.show()



    