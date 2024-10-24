import numpy as np
import matplotlib.pyplot as plt

# Define the system parameters
rho = 1.225  # air density [kg/m^3]
R = 92.4 # rotor radius [m]
Cp_opt = 0.441 # From A1 (?)
TSR_opt = 7.5 # From A1
eff = 1 # Assume an efficiency?
gearbox_ratio = 1 # Assume a gearbox ratio of 1
I_r = 0.1868745018E+09 # rotor inertia [kg m^2]
I_g = 0.1906255023E+09 # generator inertia [kg m^2]
damping_ratio = 0.7 # damping ratio of 0.7
w_n = 0.05 * 2* np.pi # Natural frequency [rad/s]

