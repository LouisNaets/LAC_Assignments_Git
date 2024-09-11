import numpy as np
import math as m
import matplotlib.pyplot as plt

aoa_design = np.deg2rad(5) #[rad]
B = 3 #[-]
TSR = 8 #[-]

r_R = np.arange(0.1, 1.01, 0.01)

Cl_design = (2*m.pi)**2/360 * aoa_design + 0.452

c_R = (16*m.pi/9) * (1 / (Cl_design * TSR**2 * B)) * 1/r_R
theta_R = 2/3 * (1/TSR * 1/r_R) - aoa_design # [rad]

# Plot c/r and theta/r as a function of r/R in subplots
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(r_R, c_R)
plt.xlabel('r/R [-]')
plt.ylabel('c/R [-]')
plt.title('Chord distribution')
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(r_R, np.rad2deg(theta_R))
plt.xlabel('r/R [-]')
plt.ylabel('Twist [deg]')
plt.title('Twist distribution')
plt.grid()
plt.show()

