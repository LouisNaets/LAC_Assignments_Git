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


# analysis settings
HAWC2S_PATH = './A4 Design Loads and AEP/Assignment/dtu_10mw_flex_minrotspd.opt'  # path to .pwr or .opt file
STATS_PATH = './A4 Design Loads and AEP/Assignment/dtu_10mw_turb_stats.hdf5'  # path to mean steady stats
SUBFOLDER = 'tca'  # which subfolder to plot: tca or tcb

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

# load the HAWC2 data from the stats file. Isolate the simulations with no tilt.
df, wsps = load_stats(STATS_PATH, subfolder=SUBFOLDER, statstype='turb')

# load/calc. the stuff we need from the HAWC2S opt/pwr file for the operational data comparisons
opt_dict = load_oper(HAWC2S_PATH)
h2s_u, h2s_pitch, h2s_rotspd, = opt_dict['ws_ms'], opt_dict['pitch_deg'], opt_dict['rotor_speed_rpm']
h2s_paero, h2s_thrust = opt_dict['power_kw'], opt_dict['thrust_kn']
h2s_aerotrq = h2s_paero / (h2s_rotspd * np.pi / 30)

# get hawc2 thrust and aerodynamic torque for theoretical calculations
h2_thrust = df.filter_channel('Thrust', CHAN_DESCS)['mean']
h2_aero_trq = df.filter_channel('GenTrq', CHAN_DESCS)['mean'] / GENEFF * 1e-3  # aerodynamic torque [kNm]

# initialize the figure and axes
fig, axs = plt.subplots(4, 3, figsize=(12, 8), clear=True)

# Set the opacity and marker size variables
dot_opacity = 0.25  # Opacity for individual points
line_opacity = 0.8  # Opacity for the mean lines
dot_size = 10       # Size for individual points
line_size = 40      # Size for mean points
max_color = 'tab:gray'
mean_color = 'tab:blue'
min_color = 'tab:orange'

# Loop over each channel and plot the steady state with the theory line
for iplot, chan_id in enumerate(chan_ids):
    
    # Isolate the channel data
    chan_df = df.filter_channel(chan_id, CHAN_DESCS)

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
    # Individual 'mean' values with lower opacity and smaller markers
    ax.scatter(h2_wind_sorted, HAWC2val_mean_sorted, color=mean_color, alpha=dot_opacity, s=dot_size, label='Mean')
    # Mean of 'mean' with higher opacity and larger markers
    ax.plot(unique_wind_speeds, mean_HAWC2val_per_wind, color=mean_color, alpha=line_opacity, markersize=line_size)#, label='Mean of means')

    # Individual 'min' values
    ax.scatter(h2_wind_sorted, HAWC2val_min_sorted, color=min_color, alpha=dot_opacity, s=dot_size, label='Min')
    # Mean of 'min'
    ax.plot(unique_wind_speeds, min_HAWC2val_per_wind, color=min_color, alpha=line_opacity, markersize=line_size)#, label='Mean of min')

    # Individual 'max' values
    ax.scatter(h2_wind_sorted, HAWC2val_max_sorted, color=max_color, alpha=dot_opacity, s=dot_size, label='Max')
    # Mean of 'max'
    ax.plot(unique_wind_speeds, max_HAWC2val_per_wind, color=max_color, alpha=line_opacity, markersize=line_size)#, label='Mean of max')

    # Formatting the plot
    ax.grid('on')
    ax.set(xlabel='Wind speed [m/s]' if iplot > 8 else None,
           ylabel=f'{chan_id} [{chan_df.units.iloc[0]}]', xlim=[4, 25])

# Add legends and format the figure
axs[0, 0].legend()
#axs[1, 2].legend()
fig.suptitle(f'Case: DTU 10MW turbine - {SUBFOLDER}')
fig.tight_layout()

plt.savefig('./A4 Design Loads and AEP/Assignment/Figures/dtu10mw_tca.svg', format='svg')
plt.savefig('./A4 Design Loads and AEP/Assignment/Figures/dtu10mw_tca.png', format='png')