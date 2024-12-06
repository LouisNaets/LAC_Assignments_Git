from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lacbox.io import ReadHAWC2
from lacbox.test import test_data_path

def plot_wind_turbine_data(data_list, omega_zeta_list, title: str, alphas=None, linestyles=None):
    """
    Plots wind turbine data for up to 4 datasets with customizable alpha (opacity) and linestyles.

    Parameters:
        data_list: List of dictionaries containing time-series data for c1, c2, c3, c4.
        omega_zeta_list: List of [omega, zeta] pairs for labeling the plots.
        title: Title of the entire figure.
        alphas: List of opacities for each dataset (default is 1.0 for all).
        linestyles: List of linestyles for each dataset (default is solid for all).
    """
    # Defaults for alpha and linestyle
    if alphas is None:
        alphas = [1.0] * len(data_list)  # Full opacity by default
    if linestyles is None:
        linestyles = ['-'] * len(data_list)  # Solid lines by default

    # Create a figure with 2x2 subplots
    fig, axs = plt.subplots(3,3, figsize=(18, 10),dpi=500)
    fig.tight_layout(pad=5.0)

    # Titles and labels for the subplots
    plot_info = [
        ("Wind Speed", "Wind Speed [m/s]", None),
        ("Pitch Angle", "Pitch [Deg]", None),
        ("Rotational Speed", "$\Omega$ [Rad/sec]", None),
        ("Aer. Torque", "Torque [kNm]", None),
        ("Aer. Thrust", "Thrust [kN]", None),
        ("Electrical Power", "Power [W]", None),
        ("TBSS", "Moment [kNm]", None),
        ("TBFA", "Moment [kNm]", None),
        #("IPBRM", "Moment [kNm]", None),
        #("EdgBRM", "Moment [kNm]", None),
        ("ShftTrs", "Moment [kNm]", None)
    ]

    # Data keys for each plot
    keys = ["wind_speed", "pitch", "rotational_speed", "torque", "thrust", "elec_power", "tbss", "tbfa", "shfttrs"]

    # Loop over subplots and datasets
    for idx, ax in enumerate(axs.flat):
        for i, data in enumerate(data_list):
            ax.plot(
                data["time"], 
                data[keys[idx]], 
                label=f'$ω_Ω$={omega_zeta_list[i][0]} Hz, $ζ_Ω$={omega_zeta_list[i][1]}',
                alpha=alphas[i], 
                linestyle=linestyles[i]
            )

        # Set subplot titles, labels, and limits
        ax.set_title(plot_info[idx][0])
        ax.set_ylabel(plot_info[idx][1])
        if plot_info[idx][2]:  # Set y-limits if specified
            ax.set_ylim(*plot_info[idx][2])
        ax.grid()

        # Add x-label for the bottom row
        if idx >= 2:
            ax.set_xlabel('Time [s]')

    # Add legend to the first plot only
    axs[0, 0].legend(loc='best')

    # Set the figure title
    #fig.suptitle(title, fontsize=16, y=0.93)
    plt.savefig(f'A5 Individual design/Figures 3/Controller_settings{title}.svg', format='svg')
    plt.savefig(f'A5 Individual design/Figures 3/Controller_settings{title}.png', format='png')

# Helper functions for data extraction
def extract_values_for_part_3_omit100(df, start=100, stop=None):
    # Extract columns
    time = df.data[:, 0]
    pitch = df.data[:, 3]
    rotational_speed = df.data[:, 9]
    elec_power = df.data[:, 104]
    wind_speed = df.data[:, 16]
    tbss = df.data[:, 20-1]
    tbfa = df.data[:, 19-1]
    thrust = df.data[:, 13-1]
    torque = df.data[:, 11-1]
    ipbrm = df.data[:, 29-1]
    edgbrm = df.data[:, 38-1]
    shfttrs = df.data[:, 27-1]

    # Manually define start and stop times for filtering
    start = start # Adjust this value as needed
    stop = stop   # Adjust this value as needed, or set to None for no upper limit
    
    # Apply filtering based on start and stop times
    if stop is not None:
        mask = (time >= start) & (time <= stop)
    else:
        mask = time >= start
    
    return {
        'pitch': pitch[mask],
        'rotational_speed': rotational_speed[mask],
        'elec_power': elec_power[mask],
        'wind_speed': wind_speed[mask],
        'time': time[mask],
        'tbfa': tbfa[mask],
        'tbss': tbss[mask],
        'thrust': thrust[mask],
        'torque': torque[mask],
        'ipbrm': ipbrm[mask],
        'edgbrm': edgbrm[mask],
        'shfttrs': shfttrs[mask]
    }

# Example usage
df1 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_1.hdf5")
df2 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_2.hdf5")
df3 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_3.hdf5")
df4 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_4.hdf5")
df5 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_5.hdf5")
df6 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_6.hdf5")
df7 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_7.hdf5")

omega_zeta_list = [[0.03, 0.7], [0.03, 0.8], [0.03, 0.9], [0.04, 0.8], [0.04, 0.9], [0.05, 0.8], [0.05, 0.9]]
data_list_max_wind = [
    extract_values_for_part_3_omit100(df1, 960, 1000),
    extract_values_for_part_3_omit100(df2, 960, 1000),
    extract_values_for_part_3_omit100(df3, 960, 1000),
    extract_values_for_part_3_omit100(df4, 960, 1000),
    extract_values_for_part_3_omit100(df5, 960, 1000),
    extract_values_for_part_3_omit100(df6, 960, 1000),
    extract_values_for_part_3_omit100(df7, 960, 1000)
]

data_list_max_wind_down = [
    extract_values_for_part_3_omit100(df1, 1001, 1041),
    extract_values_for_part_3_omit100(df2, 1001, 1041),
    extract_values_for_part_3_omit100(df3, 1001, 1041),
    extract_values_for_part_3_omit100(df4, 1001, 1041),
    extract_values_for_part_3_omit100(df5, 1001, 1041),
    extract_values_for_part_3_omit100(df6, 1001, 1041),
    extract_values_for_part_3_omit100(df7, 1001, 1041)
]

data_list_v_rated = [
    extract_values_for_part_3_omit100(df1, 386, 426),
    extract_values_for_part_3_omit100(df2, 386, 426),
    extract_values_for_part_3_omit100(df3, 386, 426),
    extract_values_for_part_3_omit100(df4, 386, 426),
    extract_values_for_part_3_omit100(df5, 386, 426),
    extract_values_for_part_3_omit100(df6, 386, 426),
    extract_values_for_part_3_omit100(df7, 386, 426)
]

data_list = [
    extract_values_for_part_3_omit100(df1),
    extract_values_for_part_3_omit100(df2),
    extract_values_for_part_3_omit100(df3),
    extract_values_for_part_3_omit100(df4),
    extract_values_for_part_3_omit100(df5),
    extract_values_for_part_3_omit100(df6),
    extract_values_for_part_3_omit100(df7)
]


#alphas = [1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1]  # Varying opacities for each dataframe
alphas = None
linestyles = None  # You can update linestyles if needed

plot_wind_turbine_data(data_list_max_wind, omega_zeta_list, "_25ms", alphas, linestyles)
plot_wind_turbine_data(data_list_max_wind_down, omega_zeta_list, "_25ms_downstep", alphas, linestyles)
plot_wind_turbine_data(data_list_v_rated, omega_zeta_list, "_11ms", alphas, linestyles)
plot_wind_turbine_data(data_list, omega_zeta_list, "", alphas, linestyles)
#plt.show()