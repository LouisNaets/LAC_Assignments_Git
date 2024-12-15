# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from lacbox.io import load_stats

# Load data for each turbine design
def load_design_data(stats_path, chan_ids, chan_descs):
    # Load stats and extract data
    df, _ = load_stats(stats_path, statstype='turb')
    data = {}
    for chan_id in chan_ids:
        chan_df = df.filter_channel(chan_id, chan_descs)
        wind_speeds = np.array(chan_df['wsp'])
        max_values = np.array(chan_df['max'])
        min_values = np.array(chan_df['min'])
        mean_values = np.array(chan_df['mean'])

        # Organize data by wind speed
        unique_ws = np.unique(wind_speeds)
        max_per_ws = [np.mean(max_values[wind_speeds == wsp]) for wsp in unique_ws]
        min_per_ws = [np.mean(min_values[wind_speeds == wsp]) for wsp in unique_ws]
        mean_per_ws = [np.mean(mean_values[wind_speeds == wsp]) for wsp in unique_ws]

        data[chan_id] = {
            'wind_speeds': unique_ws,
            'max': np.array(max_per_ws),
            'min': np.array(min_per_ws),
            'mean': np.array(mean_per_ws)
        }
    return data

# Paths to data files
chan_ids = ['TbFA', 'TbSS', 'YbTilt', 'YbRoll', 'ShftTrs', 'OoPBRM', 'IPBRM', 'EdgBRM', 'FlpBRM']
chan_descs = {
    'TbFA': 'momentmx mbdy:tower nodenr:   1',
    'TbSS': 'momentmy mbdy:tower nodenr:   1',
    'YbTilt': 'momentmx mbdy:tower nodenr:  11',
    'YbRoll': 'momentmy mbdy:tower nodenr:  11',
    'ShftTrs': 'momentmz mbdy:shaft nodenr:   4',
    'OoPBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: hub1',
    'IPBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: hub1',
    'EdgBRM': 'momentmy mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped',
    'FlpBRM': 'momentmx mbdy:blade1 nodenr:   1 coo: blade1  blade1 root flped'
}

dtu_path = './A4 Design Loads and AEP/stats_files/dtu_10mw_turb_stats.hdf5'
group7_path = './A4 Design Loads and AEP/stats_files/group7_turbB_stats_our.csv'
redesign_path = './A4 Design Loads and AEP/stats_files/group7_turbB_stats.csv'

# Load data for all designs
dtu_data = load_design_data(dtu_path, chan_ids, chan_descs)
group7_data = load_design_data(group7_path, chan_ids, chan_descs)
redesign_data = load_design_data(redesign_path, chan_ids, chan_descs)

# Define a nice color scheme
colors = ['indigo', 'mediumslateblue', 'skyblue']  # Matplotlib's viridis colormap
color_dtu = colors[0]  # First color in viridis
color_group7 = colors[1]  # Second color in viridis
color_redesign = colors[2]  # Third color in viridis

# Plot comparisons
fig, axs = plt.subplots(3, 3, figsize=(18, 15))
axs = axs.flatten()

for i, chan_id in enumerate(chan_ids):
    ax = axs[i]

    # DTU
    ax.plot(dtu_data[chan_id]['wind_speeds'], dtu_data[chan_id]['max'], marker='o', linestyle='-', color=color_dtu, label='DTU Max')
    ax.plot(dtu_data[chan_id]['wind_speeds'], dtu_data[chan_id]['min'], marker='o', linestyle='--', color=color_dtu, label='DTU Min')

    # Group7
    ax.plot(group7_data[chan_id]['wind_speeds'], group7_data[chan_id]['max'], marker='s', linestyle='-', color=color_group7, label='Group7 Max')
    ax.plot(group7_data[chan_id]['wind_speeds'], group7_data[chan_id]['min'], marker='s', linestyle='--', color=color_group7, label='Group7 Min')

    # Redesigned
    ax.plot(redesign_data[chan_id]['wind_speeds'], redesign_data[chan_id]['max'], marker='^', linestyle='-', color=color_redesign, label='Redesign Max')
    ax.plot(redesign_data[chan_id]['wind_speeds'], redesign_data[chan_id]['min'], marker='^', linestyle='--', color=color_redesign, label='Redesign Min')

    ax.set_title(chan_id, fontsize=16)
    ax.set_xlabel('Wind Speed [m/s]', fontsize=12)
    ax.set_ylabel('Load [kNm]', fontsize=12)
    ax.grid(True)
    ax.legend(fontsize=10)

# Adjust layout and save the figure
fig.tight_layout()
plt.savefig('./A4 Design Loads and AEP/figures/comparison_design_loads_updated.png', format='png')
plt.savefig('./A4 Design Loads and AEP/figures/comparison_design_loads_updated.svg', format='svg')
plt.show()
