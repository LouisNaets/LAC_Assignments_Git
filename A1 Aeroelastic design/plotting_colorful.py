import numpy as np
import matplotlib.pyplot as plt



#Variables
Cl = 1
TSR  = 10
B = 10
r = 1
R = 100
alpha = 0

#Equation A
chord_per_radius = (16*np.pi/9)*(Cl*TSR**2*B)*(1/(r/R))

#Equation B
twist_per_radius = (2/3)*(1/TSR*1/(r/R))-alpha


