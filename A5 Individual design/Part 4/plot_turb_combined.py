# -*- coding: utf-8 -*-
"""Compare mean value of steady simulations in HAWC2 (blue dots) to "theory" (lines).

For the operational parameters (i.e., pitch, rotor speed, etc.), the "theoretical" values
are the corresponding values in a .pwr file.

For the load channels of interest, the theoretical lines are calculated according to the
theoretical equations we derived for each load channel, as a function of thrust/torque/
gravity moment.

YOUR TASK! Add the lines that calculate the theory, as prompted by the slides.
"""
from lacbox.io import load_stats, load_oper, load_st
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import statistics as stats

plt.rcParams.update({'axes.labelsize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 8, 'axes.titlesize': 15})

# analysis settings
HAWC2S_PATH = './hawc_files/our_design/data/group7_3B_design_flex.opt'  # path to .pwr or .opt file
STATS_PATH = './A4 Design Loads and AEP/Assignment/group7_turbB_stats_ts.csv'  # path to mean steady stats
STATS_PATH_DTU = './A4 Design Loads and AEP/Assignment/dtu_10mw_turb_stats.hdf5'
SUBFOLDER_our = 'tcb'
SUBFOLDER_DTU = 'tca'  # which subfolder to plot: tca or tcb

# turbine constants
GENEFF = 0.94  # generator/gearbox efficienty [%]
FG_TIMES_DY = 6250  # yaw-bearing pitch moment due to gravity [kNm]
if 'notilt' in SUBFOLDER_DTU:
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
              'OoPHub': 'momentmx mbdy:hub1 nodenr:   1 coo: hub1',
              'IPHub': 'momentmy mbdy:hub1 nodenr:   1 coo: hub1',
              'FlpBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped',
              'EdgBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped',
              'TowerClearance': 'min. distance bladetips tower'
              }


#'min. distance bladetips tower'
#'TowerClerance': 'DLL :  5 inpvec :   1  min. distance bladetips tower [m]' tower clearance

#'EdgBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped' edgewise blade moment

#'FlpBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped' flapwise blade moment



# what channels we want to plot
chan_ids = ['BldPit', 'RotSpd', 'Thrust', 'GenTrq', 'ElPow', 'TbFA', 'TbSS',
            'YbTilt', 'YbRoll', 'ShftTrs', 'OoPBRM', 'IPBRM', 'EdgBRM', 'FlpBRM', 'TowerClearance']

turb_ids = ['path', 'filename', 'subfolder', 'ichan', 'names', 'units', 'desc',
            'mean', 'max', 'min', 'std', '1%', '50%', '99%', 'del3', 'del4', 'del5',
            'del8', 'del10', 'del12', 'wsp']

# load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
df, wsps = load_stats(STATS_PATH, statstype='turb')
df_DTU, wsps_DTU = load_stats(STATS_PATH_DTU, subfolder=SUBFOLDER_DTU, statstype='turb')

dfs = [[df, wsps, 0, 'Group 7'],[df_DTU, wsps_DTU, 1, 'DTU 10MW']]

# initialize the figure and axes
fig, axs = plt.subplots(5, 3, figsize=(12, 13), clear=True, dpi=500)

# Set the opacity and marker size variables
dot_opacity = 0.25  # Opacity for individual points
dot_size = 7       # Size for individual points
line_opacity = 1  # Opacity for the mean lines
line_size = 40      # Size for mean points
color_mean = ['tab:blue','tab:red']
color_outer_bounds = ['cornflowerblue','lightcoral']

for df, wsps, i, label_name in dfs:
    # Loop over each channel and plot the steady state with the theory line
    for iplot, chan_id in enumerate(chan_ids):
        
        # Isolate the channel data
        chan_df = df.filter_channel(chan_id, CHAN_DESCS)
        chan_df_DTU = df_DTU.filter_channel(chan_id, CHAN_DESCS)

        # Extract HAWC2 wind and the stats ('mean', 'min', 'max') for the channel
        h2_wind = np.array(chan_df['wsp'])
        HAWC2val_mean = np.array(chan_df['mean'])
        HAWC2val_min = np.array(chan_df['min'])
        HAWC2val_max = np.array(chan_df['max'])

        # Sort HAWC2 values by increasing wind speed
        i_h2 = np.argsort(h2_wind)
        h2_wind_sorted = h2_wind[i_h2]
        HAWC2val_mean_sorted = HAWC2val_mean[i_h2]
        HAWC2val_min_sorted = HAWC2val_min[i_h2]
        HAWC2val_max_sorted = HAWC2val_max[i_h2]

        # Calculate the mean for each unique wind speed for 'mean', 'min', and 'max'
        unique_wind_speeds = np.unique(h2_wind_sorted)
        mean_HAWC2val_per_wind = [np.mean(HAWC2val_mean_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
        min_HAWC2val_per_wind = [np.mean(HAWC2val_min_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]
        max_HAWC2val_per_wind = [np.mean(HAWC2val_max_sorted[h2_wind_sorted == wsp]) for wsp in unique_wind_speeds]

        # Plot the results
        ax = axs.flatten()[iplot]

        if chan_id == 'TowerClearance':
            # Individual 'min' values
            ax.scatter(h2_wind_sorted, HAWC2val_min_sorted, color=color_outer_bounds[i], alpha=dot_opacity, s=dot_size)
            # Mean of 'min'
            ax.plot(unique_wind_speeds, min_HAWC2val_per_wind, color=color_outer_bounds[i], alpha=line_opacity, markersize=line_size, label=f'$Min/Max_{{{label_name}}}$')#, label='Mean of min')

            # Formatting the plot
            ax.grid('on')
            ax.set(xlabel='Wind speed [m/s]' if iplot > 11 else None,
                ylabel=f'{chan_id} [m]', xlim=[4, 25])

        else:
            # Individual 'mean' values with lower opacity and smaller markers
            ax.scatter(h2_wind_sorted, HAWC2val_mean_sorted, color=color_mean[i], alpha=dot_opacity, s=dot_size)
            # Mean of 'mean' with higher opacity and larger markers
            ax.plot(unique_wind_speeds, mean_HAWC2val_per_wind, color=color_mean[i], alpha=line_opacity, markersize=line_size, label=f'$μ_{{{label_name}}}$')#, label='Mean of means')

            # Individual 'min' values
            ax.scatter(h2_wind_sorted, HAWC2val_min_sorted, color=color_outer_bounds[i], alpha=dot_opacity, s=dot_size)
            # Mean of 'min'
            ax.plot(unique_wind_speeds, min_HAWC2val_per_wind, color=color_outer_bounds[i], alpha=line_opacity, markersize=line_size, label=f'$Min/Max_{{{label_name}}}$')#, label='Mean of min')

            # Individual 'max' values
            ax.scatter(h2_wind_sorted, HAWC2val_max_sorted, color=color_outer_bounds[i], alpha=dot_opacity, s=dot_size)
            # Mean of 'max'
            ax.plot(unique_wind_speeds, max_HAWC2val_per_wind, color=color_outer_bounds[i], alpha=line_opacity, markersize=line_size)#, label='Mean of max')

            # Formatting the plot
            ax.grid('on')
            ax.set(xlabel='Wind speed [m/s]' if iplot > 11 else None,
                ylabel=f'{chan_id} [{chan_df.units.iloc[0]}]', xlim=[4, 25])
            
            if chan_id == 'GenTrq':
                ax.set(ylabel=f'{chan_id} [Nm]')
            elif chan_id == 'ElPow':
                ax.set(ylabel=f'{chan_id} [W]')

# Add legends and format the figure
axs[0, 0].legend()
#fig.suptitle(f'Case: Group 7 design - {SUBFOLDER_our}')
fig.tight_layout()

plt.savefig('./A4 Design Loads and AEP/Assignment/Figures/combined_turb_stats.svg', format='svg')
plt.savefig('./A4 Design Loads and AEP/Assignment/Figures/combined_turb_stats.png', format='png')
