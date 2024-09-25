# -*- coding: utf-8 -*-
"""Eigenanalysis for Ronnie.

Prints the modal components and, optionally, visualizes the mode shapes.
"""
import numpy as np
from scipy.linalg import eig
from _ronnie_utils import get_ff_matrices, print_modal_components, animate_mode_shape


# define our input variables
Omega = 0.8  # rotational frequency [rad/s]
animate_modes = True  # visualize the modes?
save_gif = False  # save the animation as a .gif? (requires matplotlib > 3.2.1, e.g. 3.3.1)

# =======================================================================================
# you shouldn't need to change anything below this line :)

# keyword arguments
kwargs = dict(dofs=[1, 1, 1, 1],  # flags to enable degrees of freedom [Bld1, Bld2, Bld3, x]
              M=50_000,  # equivalent rotor/nacelle/tower mass [N/m]
              kx=200_000,  # equivalent tower stiffness [N/m]
              m=500,  # equivalent blade mass [kg]
              l=30,  # equivalent blade length [m]
              omega0=2*np.pi*0.8)  # blade natural frequency [rad/s]

# get the system matrices in the fixed frame
M_ff, C_ff, K_ff = get_ff_matrices(Omega, **kwargs)

# assemble the state-space matrices
n = M_ff.shape[0]  # number of dofs
I, O = np.eye(n), np.zeros((n, n))
A_ff = np.block([[O, I],
                 [-K_ff, -C_ff]])
B_ff = np.block([[I, O],
                 [O, M_ff]])

# solve the generalized eigenvalue problem (GEVP)
[L, P] = eig(A_ff, B_ff)

# isolate the states that correspond to displacements
L = L[::2] # Every other value is removed (because of complex conjugate pairs only every other is unique)
P = P[:n, ::2] # In lines 0-n, every other value is removed

# sort in increasing order
argsort = np.argsort(np.abs(L))
L = L[argsort]
P = P[:, argsort]

# calculate natural frequencies and damping
omegas = np.abs(L)  # rad/s
freqs = omegas / 2 / np.pi  # Hz
zetas = -np.cos(np.angle(L)) * 100  # % critical

# extract the modal components in the non-rotating frame
a0, a1, b1 = P[:3, :]  # extract the blade components from eigenvectors
x = P[3, :]  # center mass displacement

# calculate the symmetric, forwards whirling and backwards wirling magnitudes
re, im = np.real, np.imag  # assign aliases to functions for convenience
S_mag = np.abs(a0)  # symmetric is just a0
BW_mag = 0.5*np.sqrt((im(a1)-re(b1))**2 + (re(a1)+im(b1))**2)  # BW, eqn in notes
FW_mag = 0.5*np.sqrt((im(a1)+re(b1))**2 + (re(a1)-im(b1))**2)  # FW, eqn in notes

# print a selection of values to screen
print_modal_components(Omega, freqs, zetas, x, a0, a1, b1, S_mag, BW_mag, FW_mag)

# visualize the mode shape, if requested
if animate_modes:
    anim = animate_mode_shape(Omega, freqs, x, a0, a1, b1)
    if save_gif:
        if Omega == 0:
            anim.save('ronnie_standstill.gif', dpi=150)
        else:
            anim.save('ronnie_rotating.gif', dpi=150)

