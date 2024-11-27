"""Visualize power-spectral density (PSD) of different load channels.

To reduce noise, PSDs are averaged across different seeds with same wind speed.


"""
from pathlib import Path

from lacbox.io import ReadHAWC2
import matplotlib.pyplot as plt
import numpy as np


def add_vert_line(ax, xval, c=None, label=None, linestyle=None, **kwargs):
    """Add vertical line to an axis at given xvalue."""
    ylim = ax.get_ylim()  # current ylim
    ax.plot([xval, xval], ylim, c=c, label=label, zorder=-1, linestyle=linestyle, **kwargs)  # plot the line
    ax.set_ylim(ylim)  # re-set y limits


# variables you can mess with
plotlog = True  # whether y axis should be logarithmic or not
plot_shftrsn = False  # whether to plot shaft torsion or not
savefig = False  # whether to save the figures to file

# define constants (do not change)
F_TOWER_MODE = 0.256  # 1st tower mode, FA and SS [Hz]
FFLAP_FIXFREE_BLADE = 0.60  # flapwise frequency of cantilever blade [Hz]
FEDGE_FIXFREE_BLADE = 0.93  # edgewise frequency of cantilever blade [Hz]
F_DRIVETRAIN = 1.82  # frequency of drivetrain mode [Hz]
CHANS = ['tbfa', 'tbss', 'oopbrm', 'ipbrm', 'shaft']  # channels we're interested in
XLIM = [0, 2.0]  # x-plot limits in Hz

# loop over wind speeds
for wsp in [6, 11, 24]:

    # ======================= CALCULATE PSDS =======================

    # get list of files corresponding to the wind speed using pattern X.X
    files = [f for f in Path('./code snippets/plot_psd/res_turb').glob('*') if f'{wsp:.1f}' in f.as_posix()]
    nfiles = len(files)

    # loop over the seeds
    for ifile, f in enumerate(files):

        # load the time-series file and get indices of relevant channels
        h2res = ReadHAWC2(f)
        names, units, desc = h2res.chaninfo
        idcs = {}
        idcs['tbfa'], idcs['tbss'] = np.where(['tower base' in d for d in desc])[0][:2]
        idcs['oopbrm'], idcs['ipbrm'] = np.where(['blade 1 root ipop' in d for d in desc])[0][:2]
        idcs['shaft'] = np.where(['shaft  main bearing' in d for d in desc])[0][-1]
        idx_omega = np.where(['Omega' in n for n in names])[0][0]
        nt = h2res.data.shape[0]
        dt = h2res.t[1] - h2res.t[0]

        # initialize arrays if first time series
        if ifile == 0:
            rot_speeds = np.empty((nfiles, nt))  # nseeds x nt array of rotor speeds
            fft_freq = np.fft.rfftfreq(nt, dt)  # array of frequencies for PSD
            nfreq = fft_freq.size  # number of frequencies in FFT vectory
            PSDs = np.empty((nfiles, nfreq, len(CHANS)))  # 3D array of FFTs, nfiles x nfreq x nchans

        # save rotor speed time series
        rot_speeds[ifile] = h2res.data[:, idx_omega]

        # calculate PSD of time series, save in array
        for ic, chan in enumerate(CHANS):
            x = h2res.data[:, idcs[chan]]  # time series
            X = np.fft.rfft(x) / nt  # FFT
            PSDs[ifile, :, ic] = 2 * nt * dt * np.abs(X)**2  # one-sided PSD: 2 abs(X)^2 / df
        PSDs[:, 0, :] = np.nan  # remove 0-Hz component (corresponds to mean)

    # sum PSDs over realizations
    PSD_avg = PSDs.sum(axis=0)  # nfreq x nchans

    # calculate rotor speed
    Omega = rot_speeds.mean() / 2 / np.pi  # in Hz

    # ======================= PLOT THINGS =======================

    # Create figure for PSDs
    fig, axs = plt.subplots(2 + plot_shftrsn, 1, figsize=(9, [6, 9][plot_shftrsn]))

    # ==== AX0: tower-base fore-aft and side-side ====
    ax = axs[0]
    ax.plot(fft_freq, PSD_avg[:, 0], label='TbFA')
    ax.plot(fft_freq, PSD_avg[:, 1], label='TbSS')
    if plotlog:
        ax.set_yscale('log')
    ax.set(xlim=XLIM)

    # add vertical lines for identified peaks (natural or forced frequencies)
    add_vert_line(ax, F_TOWER_MODE, c='g', label='1st tower mode')
    add_vert_line(ax, FEDGE_FIXFREE_BLADE-Omega, c='b', linestyle=':', label='1st edgewise (BW)')
    add_vert_line(ax, FEDGE_FIXFREE_BLADE+Omega, c='b', linestyle='--', label='1st edgewise (FW)')
    add_vert_line(ax, F_DRIVETRAIN, c='b', label='1st edgewise (SYM)')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE-Omega, c='y', linestyle=':', label='1st flapwise (BW)')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE+Omega, c='y', linestyle='--', label='1st flapwise (FW)')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE, c='y', label='1st flapwise (SYM)')
    add_vert_line(ax, 1*Omega, c='r', linestyle=':', label='1P-harmonics')
    add_vert_line(ax, 3*Omega, c='r', linestyle='--', label='3P-harmonics')
    add_vert_line(ax, 6*Omega, c='r', label='6P-harmonics')

    ax.grid()
    ax.legend()

    # ==== AX1: blade out-of-plane and in-plane blade-root moments ====
    ax = axs[1]
    ax.plot(fft_freq, PSD_avg[:, 2], label='OoPBRM')
    ax.plot(fft_freq, PSD_avg[:, 3], label='IPBRM')
    if plotlog:
        ax.set_yscale('log')
    ax.set(xlim=XLIM)

    # add vertical lines for identified peaks (natural or forced frequencies)
    # ... !!!TODO!!! Add lines below using constants defined at top of script
    # (add lines here)
    add_vert_line(ax, F_TOWER_MODE, c='g', label='Tower')
    add_vert_line(ax, F_TOWER_MODE-Omega, c='g', linestyle='--')#, label='Tower(-)')
    add_vert_line(ax, F_TOWER_MODE+Omega, c='g', linestyle=':')#, label='Tower(+)')

    add_vert_line(ax, FFLAP_FIXFREE_BLADE, c='y', label='Flapwise FW/BW/SYM')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE-Omega, c='y', linestyle='--')#, label='Flapwise FW/BW/SYM(-)')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE-2*Omega, c='y', linestyle='--')#, label='Flapwise FW/BW/SYM(--)')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE+Omega, c='y', linestyle=':')#, label='Flapwise FW/BW/SYM(+)')
    add_vert_line(ax, FFLAP_FIXFREE_BLADE+2*Omega, c='y', linestyle=':')#, label='Flapwise FW/BW/SYM(++)')

    add_vert_line(ax, FEDGE_FIXFREE_BLADE, c='r', label='Edgewise FW/BW')
    add_vert_line(ax, FEDGE_FIXFREE_BLADE-Omega, c='r', linestyle='--')#, label='Edgewise FW/BW(-)')
    add_vert_line(ax, FEDGE_FIXFREE_BLADE-2*Omega, c='r', linestyle='--')#, label='Edgewise FW/BW(--)')
    add_vert_line(ax, FEDGE_FIXFREE_BLADE+Omega, c='r', linestyle=':')#, label='Edgewise FW/BW(+)')
    add_vert_line(ax, FEDGE_FIXFREE_BLADE+2*Omega, c='r', linestyle=':')#, label='Edgewise FW/BW(++)')

    add_vert_line(ax, F_DRIVETRAIN, c='b', label='Edgewise SYM')
    add_vert_line(ax, F_DRIVETRAIN, c='b', linestyle='--')#, label='Edgewise SYM(-)')
    add_vert_line(ax, F_DRIVETRAIN, c='b', linestyle=':')#, label='Edgewise SYM(+)')

    add_vert_line(ax, 1*Omega, c='purple', linestyle=':', label='1P-harmonics')
    add_vert_line(ax, 3*Omega, c='purple', linestyle='--', label='3P-harmonics')
    add_vert_line(ax, 6*Omega, c='purple', label='6P-harmonics')  

    ax.grid()
    ax.legend(ncols=2)

    # ==== (if requested) AX2: shaft torsion ====
    if plot_shftrsn:
        ax = axs[2]
        ax.plot(fft_freq, PSD_avg[:, 4], label='Shaft')
        if plotlog:
            ax.set_yscale('log')
        ax.set(xlim=XLIM)

    # add vertical lines for identified peaks (natural or forced frequencies)
    # ... Optional: Add lines below using constants defined at top of script
    # (add lines here)

        ax.grid()
        ax.legend(ncols=2)

    fig.suptitle(f'Wind speed: {wsp} m/s   1P: {Omega:.2f} Hz')
    fig.tight_layout()
    if savefig:
        fig.savefig(f'psd_{wsp:.1f}.png', dpi=150)

plt.show()

print()