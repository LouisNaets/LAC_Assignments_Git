from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lacbox.io import ReadHAWC2
from lacbox.test import test_data_path

plt.rcParams.update({'axes.labelsize': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 10, 'axes.titlesize': 15})

def plot_wind_turbine_data(c1,c2, title:str, omega1:str, omega2:str, zeta1:str, i:int):
    # Create a figure and two sets of 2x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))

    # Adjust the layout to prevent overlapping
    #fig.subplots_adjust(hspace=0.2, wspace=0.1)
    fig.tight_layout(pad=5.0)
    #fig.suptitle(f"{title}", fontsize=18, y=0.93)

    # C1-3 plots (first set)
    axs[0, 0].plot(c1["time"], c1["wind_speed"], label=f'C{i+1}: $ω_Ω$=0.05 Hz & $ζ_Ω=0.7 $')
    axs[0, 0].plot(c2["time"], c2["wind_speed"], alpha=0.65, label=f'C{i+7}: $ω_Ω$={omega1} Hz & $ζ_Ω={zeta1} $')
    axs[0, 0].set_title('Wind Speed', fontsize=18)
    axs[0, 0].set_ylabel('Wind Speed [m/s]', fontsize=18)
    axs[0, 0].legend(loc='best', fontsize=14)
    axs[0, 0].tick_params(axis='both', which='major', labelsize=16)
    axs[0, 0].tick_params(labelbottom=False)  # Remove x-axis tick marks
    axs[0, 0].grid()

    axs[0, 1].plot(c1["time"], c1["pitch"])
    axs[0, 1].plot(c2["time"], c2["pitch"], alpha=0.65)
    axs[0, 1].set_title('Pitch Angle', fontsize=18)
    axs[0, 1].set_ylabel('Pitch [Deg]', fontsize=18)
    axs[0, 1].yaxis.set_label_position("right")
    axs[0, 1].yaxis.tick_right()
    axs[0, 1].set_ylim(-1, 25)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=16)
    axs[0, 1].tick_params(labelbottom=False)  # Remove x-axis tick marks
    axs[0, 1].grid()

    axs[1, 0].plot(c1["time"], c1["rotational_speed"])
    axs[1, 0].plot(c2["time"], c2["rotational_speed"], alpha=0.65)
    axs[1, 0].set_title('Rotational Speed', fontsize=18)
    axs[1, 0].set_xlabel('Time [s]', fontsize=18)
    axs[1, 0].set_ylabel('$\Omega$ [Rad/sec]', fontsize=18)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=16)
    axs[1, 0].grid()

    axs[1, 1].plot(c1["time"], c1["elec_power"], label='$ω_Ω$=0.05 Hz & $ζ_Ω=0.7 $')
    axs[1, 1].plot(c2["time"], c2["elec_power"], alpha=0.65, label=f'$ω_Ω$={omega1} Hz & $ζ_Ω={zeta1} $')
    axs[1, 1].set_title('Electrical Power', fontsize=18)
    axs[1, 1].set_xlabel('Time [s]', fontsize=18)
    axs[1, 1].set_ylabel('Power [W]', fontsize=18)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=16)
    axs[1, 1].yaxis.set_label_position("right")
    axs[1, 1].yaxis.tick_right()
    axs[1, 1].grid()

    plt.savefig('A3 Controller design/Assignment/Figures/c1c7.svg', format='svg')
    plt.savefig('A3 Controller design/Assignment/Figures/c1c7.png', format='png')
    
def extract_values_for_part_3(df):

    return {
        'pitch': df.data[:,3],
        'rotational_speed': df.data[:,9],
        'elec_power': df.data[:,104],
        'wind_speed': df.data[:, 16],
        'time': df.data[:,0]
    }    

def extract_values_for_part_3_omit100(df):
    # Extract columns
    time = df.data[:, 0]
    pitch = df.data[:, 3]
    rotational_speed = df.data[:, 9]
    elec_power = df.data[:, 104]
    wind_speed = df.data[:, 16]
    
    # Filter data to start from 100 seconds
    mask = time >= 100
    return {
        'pitch': pitch[mask],
        'rotational_speed': rotational_speed[mask],
        'elec_power': elec_power[mask],
        'wind_speed': wind_speed[mask],
        'time': time[mask]
    }

df1 = ReadHAWC2("hawc_files/our_design/res/group7_3B_design_A3_part3_C1.hdf5")
df7 = ReadHAWC2("hawc_files/our_design/res/group7_3B_design_A3_part3_C7_1.hdf5")

c1 = extract_values_for_part_3_omit100(df1)
c7 = extract_values_for_part_3_omit100(df7)

plot_wind_turbine_data(c1, c7, "", '0.03', '0.03', '0.7', 0)
#plt.show()