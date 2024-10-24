# -*- coding: utf-8 -*-
"""
Part 1
Created on Thu Jan 12 14:03:08 2023

@author: wali
"""

import matplotlib
matplotlib.use('Qt5Agg')

import matplotlib as plt
import lib # importing modules from the folder
from lib.lib_main import *

Simulation_TEND = 300.0 # Simulation lenght in seconds
SimParams = SimParams_(Simulation_TEND)

# Desired outputs
Power = 10e6
Rotor_speed = 1.005

Open_Loop_Pitch = 7.824 #  pitch angle in degrees - perfect value: 7.824
Open_Gen_Torque = Power/Rotor_speed # generator reaction torque in Nm

Controller_Type = 'OL' # Controller Type, OL:Open Loop
Controller = Controller_(Controller_Type, Open_Loop_Pitch = Open_Loop_Pitch, Open_Gen_Torque  = Open_Gen_Torque)

WT_Model = 'WT0' # Model complexity,  WT0:Rotor,   WT1:Rotor+DT,  WT2:Rotor+DT+Tower fore-aft
WT = WT_(WT_Model,SimParams)

# ------------------------------------------------------------------------------------
#   Choose the wind speed profile here:
# ------------------------------------------------------------------------------------
# Here you can choose the type of wind speed you'd like to use!
# wind_profile_options = 1
# 1: for step wind speed, use WSP_Profile_Generator to produce wind steps!
# 2: for stochastic wind speed, mean wind speed: 8 m/s
# 3: for stochastic wind speed, mean wind speed: 12 m/s
# 4: for stochastic wind speed, mean wind speed: 15 m/s
# 5: for stochastic wind speed, mean wind speed: 18 m/s
# 6: for stochastic wind speed, mean wind speed: 15 m/s ,no wind shear.
# 7: EOG
# 8: for step wind speed below rated, use WSP_Profile_Generator to produce wind steps!

wind_profile_options = 1

data = simulate(SimParams, WT, Controller, wind_profile_options)

figsize = (8,8)
gen_plot(WT, SimParams, data,figsize)

plt.show()