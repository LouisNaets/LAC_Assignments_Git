# -*- coding: utf-8 -*-
"""Utility functions for the Coleman demo
"""
from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np


def animate_mode_shape(Omega, freqs, xs, a0s, a1s, b1s):
    """"""
    # variables we need for the plot
    time = np.linspace(0, 1, num=300, endpoint=False)
    dt = 1/30*1000  # milliseconds
    wk = 2*np.pi*2  # hard-code natural frequency for sake of visualiztion
    Omega = [2*np.pi*1, 0][Omega == 0]  # Omega is either hard-coded or 0
    max_defls = np.abs(a0s) + np.abs(a1s) + np.abs(b1s)  # maximum possible blade deflections
    
    # initialize figure, axes and artists
    fig, axs = plt.subplots(1, 4, num=1, clear=True, figsize=(7, 2))
    all_artists = []
    for i, ax in enumerate(axs):
        ax.axis('equal')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_title(f'Mode {i+1} ({freqs[i]:.2} Hz)', fontsize=8)
        ax.plot(0, 0, 'o', c='0.8')
        ax.axis('off')
        artists = ([ax.plot([], [], '--', lw=2, c='0.6')[0] for i in range(3)]  # undeflected blade
                   + [ax.plot([], [], lw=2)[0] for i in range(3)]  # deflected blade
                   + [ax.plot([], [], marker='o', lw=2)[0] for i in range(3)]  # blade tip
                   + [ax.plot([], [], marker='o', c='0.4', mfc='0.9', lw=2)[0]])  # COG
        all_artists = all_artists + artists
    plt.tight_layout()
    plt.figtext(0.015, 0.07, '©Rinker')

    def update(t):
        """Update the data in the artists based on time step"""
        for i, ax in enumerate(axs):  # loop over modes
            a0, a1, b1, x, max_defl = a0s[i], a1s[i], b1s[i], xs[i], max_defls[i]  # modal comps
            x0 = np.abs(x) * np.cos(wk*t + np.angle(x))  # tower oscillation
            artists = all_artists[i*10:(i+1)*10]  # isolate artists for this axis
            x_cog, y_cog = 0, 0  # initialize cog
            for i in range(3):  # loop over blades
                # isolate and calculate things we need
                origin = artists[i]  # where deflection is measured from
                blade = artists[i + 3]  # blade
                tip = artists[i + 6]  # blade tip
                Psi = Omega * t + i*2*np.pi/3  # azimuthal angle of blade origin
                theta_i = (np.abs(a0) * np.cos(wk*t + np.angle(a0))  # blade deflection
                            + np.abs(a1) * np.cos(wk*t + np.angle(a1)) * np.cos(Psi)
                            + np.abs(b1) * np.cos(wk*t + np.angle(b1)) * np.sin(Psi))
                # plot the undeflected blade position
                x = [x0, x0 + np.sin(Psi)]
                y = [0, np.cos(Psi)]
                origin.set_data(x, y)
                # calculate tip position wrt. rotor origin and update COG
                x_tip = 0.9*np.sin(Psi + theta_i)
                y_tip = 0.9*np.cos(Psi + theta_i)
                x_cog += x_tip / 3
                y_cog += y_tip / 3
                # update plotting data
                xb = [x0, x0 + x_tip]
                yb = [0, y_tip]
                blade.set_data(xb, yb)
                tip.set_data([x0 + x_tip], [y_tip])
                # color blades and tip based on their deflection
                rel_defl = (theta_i + max_defl) / (2 * max_defl)
                c = plt.get_cmap('coolwarm')(rel_defl)
                blade.set_color(c)
                blade.set_alpha(0.5)
                tip.set_color(c)
                tip.set_alpha(0.5)
            # update COG
            cog = artists[9]  # center of gravity
            scale = max_defl*20
            cog.set_data([x0 + x_cog*scale], [y_cog*scale])  # scale for easier viz.
        return all_artists
    
    # animate the plot
    anim = animation.FuncAnimation(fig, update,  # figure handle and update functions
                                   frames=time,  # array of times
                                   interval=dt,  # delay between frames [ms]
                                   blit=True)
    
    plt.show()
    
    return anim


def get_eig(M, C, K):
    """Eigenvalues for mass, stiffness and damping matrix (assumes M invertible)"""
    A = get_state_space(M, C, K)
    s, v = np.linalg.eig(A)  # standard evp
    return s, v


def get_mode_mags(v, dofs=(1, 1, 1, 1)):
    """Magnitudes of components from eigenvector(s) in non-rotating frame"""
    dofs = dofs + [0]  # hide y dof from students
    if sum(dofs[:3]) < 3:
        raise ValueError('Mode magnitude code breaks if a blade DOF is disabled!')
    re, im = np.real, np.imag  # assign handles for convenience
    ndofs = sum(dofs)
    mode_mags = np.empty((ndofs, ndofs))  # initialize array (as, bw, fw, x, y by freq)
    a0, a1, b1 = v[:3]  # extract the blade components from eigenvectors
    mode_mags[0, :] = np.abs(a0)  # symmetric is just a0
    mode_mags[1, :] = 0.5*np.sqrt((im(a1)-re(b1))**2
                                  + (re(a1)+im(b1))**2)  # BW, eqn in dynamics handbook
    mode_mags[2, :] = 0.5*np.sqrt((im(a1)+re(b1))**2
                                  + (re(a1)-im(b1))**2)  # FW
    if dofs[3]:  # if lateral motion is enabled
        mode_mags[2, :] = np.abs(v[3])  # tower is just tower
    if dofs[4]:  # if vertical motion is enabled
        mode_mags[-1, :] = np.abs(v[4])  # tower is just tower
    return mode_mags


def get_ff_matrices(Omega, dofs=(1, 1, 1, 1), M=50_000, kx=200_000,
                    ky=200_000, m=500, l=30, omega0=2*np.pi*0.8):
    """System matrices in the non-rotating frame"""
    dofs = dofs + [0]  # hide y DOF from students
    # useful intermediates
    ml = m*l
    ml2 = m*l**2
    # system matrices
    M_nr = np.array([[ml2, 0, 0, 0, 0],
                     [0, ml2, 0, ml, 0],
                     [0, 0, ml2, 0, -ml],
                     [0, 3*ml/2, 0, M+3*m, 0],
                     [0, 0, -3*ml/2, 0, M+3*m]])
    C_nr = np.array([[0, 0, 0, 0, 0],
                     [0, 0, 2*ml2*Omega, 0, 0],
                     [0, -2*ml2*Omega, 0, 0, 0],
                     [0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0]])
    K_nr = np.diag([ml2*omega0**2, ml2*(omega0**2 - Omega**2),
                    ml2*(omega0**2 - Omega**2), kx, ky])
    # single out only the degrees of freedom we want
    idxs = np.asarray(dofs, dtype=bool)
    M_nr = M_nr[idxs][:, idxs]
    C_nr = C_nr[idxs][:, idxs]
    K_nr = K_nr[idxs][:, idxs]
    return M_nr, C_nr, K_nr


def get_state_space(M, C, K):
    """State-space matrix A, assuming M is invertible"""
    Minv = np.linalg.inv(M)  # for convenience
    n = M.shape[0]  # no. dofs
    A = np.block([[np.zeros((n, n)), np.eye(n)],
                  [-Minv@K, -Minv@C]])  # state space
    return A


def plot_campbell(Omegas=np.linspace(0, 1, 21), num=None, dofs=(1, 1, 1, 1), M=50_000,
                  kx=200_000, ky=200_000, m=500, l=30, omega0=2*np.pi*0.8):
    """Plot the Campbell diagram for the system, return eigenvalues and vectors"""
    dofs = dofs + [0]  # hide y DOF from students
    # intermediate variable
    nplot, ndofs = Omegas.size, sum(dofs)
    # initialize the arrays of undamped natural frequencies and eigenvectors
    fns = np.empty((nplot, ndofs))  # array of freqs
    vs = np.empty((nplot, ndofs, ndofs), dtype=complex)
    f_BW_FW = np.empty((nplot, 2), dtype=int)  # freqs of BW/FW modes
    # loop over the rotational frequencies to make Campbell diagram
    for iOm in range(nplot):
        # get non-rotating matrices
        m_nr, c_nr, k_nr = get_nr_matrices(Omegas[iOm], dofs=dofs, M=M, kx=kx, ky=ky,
                                           m=m, l=l, omega0=omega0)  # get matrices
        # eigenvalues and eigenvectors
        s_nr, v_nr = get_eig(m_nr, c_nr, k_nr)  # calc evals
        fns[iOm, :] = np.abs(s_nr[::2]) / (2*np.pi)  # evals -> freqs
        vs[iOm, :, :] = v_nr[:ndofs, ::2]  # only disp (not vel), unique
    # make the campbell diagram
    plt.figure(num, figsize=(8, 4)); plt.clf();
    plt.plot(Omegas, fns, 'o')
    plt.ylim([0, 1])
    plt.grid('on')
    plt.xlabel('Omega [rad/s]'); plt.ylabel('Frequency [Hz]');
    plt.tight_layout()  # make the layout nice
    plt.show()
    return fns, vs


def plot_time_freq_response(t, out, num=None):
    """Plot a simulated response in the time and frequency domains"""
    T = t[-1] + t[1]  # total time to simulate
    out_fft = np.fft.rfft(out, axis=0)
    out_f = np.arange(out_fft.shape[0]) / T
    # initialize/clear the figure
    fig = plt.figure(num=num)
    plt.clf()
    # time domain subplots
    plt.subplot(2, 2, 1)  # blades
    plt.plot(t, out[:, :3])
    plt.xlim([0, T])
    plt.grid('on')
    plt.ylabel('Blade defl. [rad]')
    plt.title('Blades, Time Domain')
    plt.subplot(2, 2, 3)  # tower
    plt.plot(t, out[:, 3])
    plt.xlim([0, T])
    plt.xlabel('Time [s]')
    plt.grid('on')
    plt.ylabel('Towertop defl. [m]')
    plt.title('Tower, Time Domain')
    # frequency domain subplots
    plt.subplot(2, 2, 2)  # blades
    plt.semilogy(out_f, np.abs(out_fft[:, :3]))
    plt.xlim([0, 1.5])
    plt.grid('on')
    plt.legend(labels=['Blade 1', 'Blade 2', 'Blade 3'])
    plt.ylabel('$|X_{bld}|(f)$')
    plt.title('Blades, Frequency Domain')
    plt.subplot(2, 2, 4)  # tower
    plt.semilogy(out_f, np.abs(out_fft[:, 3]))
    plt.xlim([0, 1.5])
    plt.grid('on')
    plt.xlabel('Freq [Hz]')
    plt.ylabel('$|X_{twr}|(f)$')
    plt.title('Tower, Frequency Domain')
    # clean up the figure
    plt.tight_layout()
    plt.show()
    return fig


def print_modal_components(Omega, freqs, zetas, x, a0, a1, b1, S_mag, BW_mag, FW_mag):
    """If Omega = 0, print freqs, damping, and imag part of x, a0, a1 and b1.
    If Omega > 0, also print symmetric and whirling components.
    """
    im = np.imag  # assign aliases to functions for convenience
    print(f'\nModal components for Omega = {Omega:.2f} rad/s')
    print('----------------------------------------')
    print('          Mode1   Mode2   Mode3   Mode4')
    # frequencies and damping
    for label, val in zip(['fk (Hz)', 'z (%C) '],
                          [freqs, zetas]):
        print(f'{label}' + ''.join([f'{v:8.3f}' for v in val]))
    print()
    # direct modal components
    if np.isclose(Omega, 0):
        for label, val in zip(['Im(a0) ', 'Im(a1) ', 'Im(b1) ', 'Im(x)  '],
                              [im(a0), im(a1), im(b1), im(x)]):
            print(f'{label}' + ''.join([f'{v:8.3f}' for v in val]))
    # symmetric, backwards and forwards whirling
    if Omega > 0:
        for label, val in zip(['a0     ', 'a1     ', 'b1     ', 'x      '],
                              [a0, a1, b1, x]):
            print(f'{label}' + ''.join([f'{v:14.3f}' for v in val]))
        print()
        for label, val in zip(['abs(S) ', 'abs(BW)', 'abs(FW)', 'abs(x) ', ],
                              [S_mag, BW_mag, FW_mag, np.abs(x), ]):
            print(f'{label}' + ''.join([f'   {v:.3f}' for v in val]))
    return
