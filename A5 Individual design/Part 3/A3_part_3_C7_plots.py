from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lacbox.io import ReadHAWC2
from lacbox.test import test_data_path

def plot_wind_turbine_data(c1,c2,c3, title:str, omega1:str, omega2:str, zeta1:str):
    # Create a figure and two sets of 2x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(18, 10))

    # Adjust the layout to prevent overlapping
    fig.tight_layout(pad=5.0)
    #fig.subplots_adjust(hspace=0.2, wspace=0.1)
    #fig.suptitle(f"{title}", fontsize=16, y=0.93)
    # C1-3 plots (first set)
    axs[0, 0].plot(c1["time"], c1["wind_speed"])
    axs[0, 0].plot(c2["time"], c2["wind_speed"])
    axs[0, 0].plot(c3["time"], c3["wind_speed"])
    axs[0, 0].set_title('Wind Speed')
    # axs[0, 0].set_xlabel('Time')
    axs[0, 0].set_ylabel('Wind Speed [m/s]')
    axs[0, 0].grid()

    axs[0, 1].plot(c1["time"], c1["pitch"], label ='$ω_Ω$=0.05 Hz & $ζ_Ω=0.7 $')
    axs[0, 1].plot(c2["time"], c2["pitch"], label =f'$ω_Ω$={omega1} Hz & $ζ_Ω={zeta1} $')
    axs[0, 1].plot(c3["time"], c3["pitch"], label =f'$ω_Ω$={omega2} Hz & $ζ_Ω={zeta1} $')
    axs[0, 1].set_title('Pitch Angle')
    # axs[0, 1].set_xlabel('Time')
    axs[0, 1].legend(loc='upper left')
    axs[0, 1].set_ylabel('Pitch [Deg]')
    axs[0, 1].set_ylim(-1, 25)
    axs[0, 1].grid()


    axs[1, 0].plot(c1["time"], c1["rotational_speed"])
    axs[1, 0].plot(c2["time"], c2["rotational_speed"])
    axs[1, 0].plot(c3["time"], c3["rotational_speed"])
    axs[1, 0].set_title('Rotational Speed')
    axs[1, 0].set_xlabel('Time [s]')
    axs[1, 0].set_ylabel('$\Omega$ [Rad/sec]')
    axs[1, 0].grid()


    axs[1, 1].plot(c1["time"], c1["elec_power"])
    axs[1, 1].plot(c2["time"], c2["elec_power"])
    axs[1, 1].plot(c3["time"], c3["elec_power"])
    axs[1, 1].set_title('Electrical Power')
    axs[1, 1].set_xlabel('Time [s]')
    
    axs[1, 1].set_ylabel('Power [W]')
    axs[1,1].grid()
    
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

df1 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_1.hdf5")
df2 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_2.hdf5")
df3 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_3.hdf5")
df4 = ReadHAWC2("hawc_files/individual_design/res/individual_design_cnew_4.hdf5")

omega_zeta_list = [[0.03, 0.7], [0.03, 0.8], [0.03, 0.9], [0.05, 0.8]]
c1 = extract_values_for_part_3_omit100(df1)
c2 = extract_values_for_part_3_omit100(df2)
c3 = extract_values_for_part_3_omit100(df3)
c4 = extract_values_for_part_3_omit100(df4)

plot_wind_turbine_data(c1, c2, c3, "", '0.03', '0.02', '0.7')
plt.show()