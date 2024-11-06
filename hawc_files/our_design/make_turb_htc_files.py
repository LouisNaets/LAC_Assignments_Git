"""Create and save a series of steady-wind files for gbar.

4 different cases, each saved in its own subfolder:
    * tilt: With tilt, flexible tower/blades, aerodynamic drag.
    * notilt: No tilt, flexible tower/blades, aerodynamic drag.
    * notiltrigid: No tilt, rigid tower/blades, aerodynamic drag.
    * notiltnodragrigid: No tilt, rigid tower/blades, no aerodynamic drag.
"""
from pathlib import Path
import random

from lacbox.htc import _clean_directory
from lacbox.io import load_oper
from myteampack import MyHTC
import numpy as np


def get_initial_rotor_speed(wsp, opt_path):
    """Given a wind speed and path to opt file, find initial rotor speed.

    Args:
        wsp (int, float): Wind speed [m/s].
        opt_path (str, pathlib.Path): Path to opt file.

    Returns:
        int, float: Initial rotor speed interpolated from opt file [rad/s].
    """
    opt_dict = load_oper(opt_path)
    opt_wsps = opt_dict['ws_ms']
    opt_rpm = opt_dict['rotor_speed_rpm']
    omega_rpm = np.interp(wsp, opt_wsps, opt_rpm)
    omega0 = omega_rpm * np.pi / 30  # rpm to rad/s
    return omega0

def get_turbulence_intensity(wind_speed, turb_class):
    """Calculate turbulence intensity based on wind class and speed.
    
    Args:
        wind_speed (float): Mean wind speed [m/s].
        turb_class (str): Turbulence class ('A', 'B', or 'C').
    
    Returns:
        float: Turbulence intensity (TI).
    """
    # Define reference turbulence intensity values for each class
    Iref_dict = {'A': 0.16, 'B': 0.14, 'C': 0.12}
    
    # Get Iref for the given class
    Iref = Iref_dict.get(turb_class.upper())
    if Iref is None:
        raise ValueError("Invalid turbulence class. Choose from 'A', 'B', or 'C'.")
    
    # Compute the standard deviation of wind speed (sigma_u)
    sigma_u = Iref * (0.75 * wind_speed + 5.6)
    
    # Calculate and return turbulence intensity
    TI = sigma_u / wind_speed
    return TI


def make_single_turb(htc, wsp, turbclass, htc_dir='./htc_turb/', res_dir='./res_turb/',
                     subfolder='', opt_path=None, seed=1337, time_start=100, time_stop=700,
                     dy=190, dz=190):
    """Make a single turbulent-wind file from a master file.
    """
    nx, ny, nz = 1024, 32, 32  # hard-code turbbox size for lac course
    # define the append name based on the subfolder and the wind speed
    wsp_seed_str = ('%.1f' % (wsp)).zfill(4) + ('_%i' % seed)  # e.g., '05.0_1337'
    if subfolder:
        append = f'_turb_{subfolder}_{wsp_seed_str}'  # e.g., '_turb_tca_05.0_1337'
    else:
        append = f'_turb_{wsp_seed_str}'  # e.g.,, '_turb_05.0_1337'
    # get new filename (excl extension) from HTCFile attribute "filename"
    fname = Path(htc.filename).name.replace('.htc', append)
    # delete hawcstab2 block
    del htc.hawcstab2
    # correct initial rotor speed if opt file is given
    if opt_path is not None:
        omega0 = get_initial_rotor_speed(wsp, opt_path)
        htc._set_initial_rotor_speed(omega0)
    # set the start and stop time
    htc.set_time(start=time_start, stop=time_stop)  # simulation times
    # calculate turbulence intensity for this turbulence class and wind speed
    turb_int = get_turbulence_intensity(wsp,turbclass)
    #turbulence = turbclass.get_turbulence(wsp)
    # set parameters in wind block
    htc.wind.tint = turb_int  # set TI
    htc.wind.turb_format = 1  # set turbulence to mann
    htc.wind.tower_shadow_method = 3  # no tower shadow
    htc.wind.wsp = wsp  # mean wind speed
    htc.wind.shear_format = [3, 0.2]  # power-law shear profile
    #htc.wind.shear_format = [1, wsp]  # constant wsp profile with height
    # set parameters in mann block
    turb_filesname = [f'./turb/{fname}_turb_{c}.bin' for c in 'uvw']
    no_grid_points = (nx, ny, nz)
    box_dimension = (wsp * (time_stop - time_start), dy, dz)
    ## high_frq_compensation=1 changed from 0 to match jennis htc file.
    htc.add_mann_turbulence(L=29.4, ae23=1, Gamma=3.9,
                            seed=seed, high_frq_compensation=1,
                            filenames=turb_filesname, no_grid_points=no_grid_points,
                            box_dimension=box_dimension,
                            dont_scale=False)
    # update name and save file (reprint of _update_name_and_save() b.c. missing kwargs to set_name)
    save_dir = Path(htc_dir)  # sanitize inputs
    # set filename using HTCFile method
    htc.set_name(fname, resdir=res_dir, subfolder=subfolder, htcdir=htc_dir)
    # save the file
    htc.save((save_dir / subfolder / (fname + '.htc')).as_posix())
    return


def main():
    """Create the HTC files for different cases, adjusting settings.
    Generate HTC files for both turbulence classes A and B, with multiple random seeds per wind speed.
    """
    # Constants for this script
    del_htc_dir = True  # Delete HTC directory if it already exists?
    master_htc = './hawc_files/our_design/_master/group7_3B_design_A4.htc'
    opt_path = './hawc_files/our_design/data/group7_3B_design_flex.opt'
    wsps = range(5, 25)  # Wind speed range
    htc_dir = './htc_turb/'  # Folder to save HTC files
    res_dir = './res_turb/'  # Where HAWC2 should save results
    start_seed = 42  # Seed for reproducibility
    num_seeds_per_wsp = 6  # Number of seeds per wind speed

    # Delete the top-level directory if requested
    _clean_directory(htc_dir, del_htc_dir)
    
    # Initialize random generator
    random.seed(start_seed)

    # Loop over both turbulence classes A and B
    for turbclass in ['A', 'B']:
        subfolder = f'tc{turbclass.lower()}'
        
        # Loop over each wind speed in the specified range
        for wsp in wsps:
            
            # Generate multiple random seeds for each wind speed
            for i in range(num_seeds_per_wsp):
                sim_seed = random.randrange(int(2**16))
                
                # Create a new HTC instance from the master file
                htc = MyHTC(master_htc)
                
                # Generate and save the HTC file for the current configuration
                make_single_turb(
                    htc=htc,
                    wsp=wsp,
                    turbclass=turbclass,
                    htc_dir=htc_dir,
                    res_dir=res_dir,
                    subfolder=subfolder,
                    opt_path=opt_path,
                    seed=sim_seed
                )



# the "script" part of this file
if __name__ == '__main__':
    main()
