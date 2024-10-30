"""Create and save a series of steady-wind files for gbar.

4 different cases, each saved in its own subfolder:
    * tilt: With tilt, flexible tower/blades, aerodynamic drag.
    * notilt: No tilt, flexible tower/blades, aerodynamic drag.
    * notiltrigid: No tilt, rigid tower/blades, aerodynamic drag.
    * notiltnodragrigid: No tilt, rigid tower/blades, no aerodynamic drag.
"""
from pathlib import Path

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


# this should have been a method in MyHTC! my bad...
def make_single_steady(htc, wsp, htc_dir='./htc_steady/', res_dir='./res_steady/', subfolder='',
                       opt_path=None, tilt=None, rigid=False, withdrag=True, time_start=200,
                       time_stop=400):
    """Make a singlesteady-wind file.

    Args:
        htc (MyHTC): MyHTC object, instantiated from the master HTC file.
        wsp (int, float): Wind speed [m/s].
        htc_dir (str or pathlib.Path, optional): Top-level directory where htc file
            should be saved. Defaults to './htc_steady/'.
        res_dir (str or pathlib.Path, optional): Top-level directory where res files
            should be saved. Defaults to './res_steady/'.
        subfolder (str or pathlib.Path, optional): Subfolder for all htc, res, and log
            files. Defaults to '' (no subfolder).
        opt_path (str or pathlib.Path, optional): Option to give the path to the opt
            file to set the initial rotor speed. Defaults to None (initial rotor speed
            not changed).
        tilt (int or float, optional): Tilt angle. Defaults to None, i.e., keep current
            tilt angle.
        rigid (bool, optional): Rigid tower and blades. Defaults to False.
        withdrag (bool, optional): Include aerodynamic drag? Defaults to True.
        time_start (int or float, optional): Time to start recording. Defaults to 200.
        time_stop (int or float, optional): Time to stop recording. Defaults to 400.

    Raises:
        ValueError: _description_
    """
    # delete hawcstab2 block
    del htc.hawcstab2
    # correct initial rotor speed if opt file is given
    if opt_path is not None:
        omega0 = get_initial_rotor_speed(wsp, opt_path)
        htc._set_initial_rotor_speed(omega0)
    # set the start and stop time
    htc.set_time(start=time_start, stop=time_stop)  # simulation times
    # update tilt if a number is passed in
    if type(tilt) in [int, float]:
        shaft_ori = htc.new_htc_structure.orientation.relative__2
        shaft_ori.mbdy2_eulerang__2 = [tilt, 0, 0]
    elif tilt is not None:
        raise ValueError('Keyword argument "tilt" must be None, int or float!')
    # rigid tower/blades if requested
    if rigid:
        htc.new_htc_structure.main_body.timoschenko_input.set = [1, 2]
        htc.new_htc_structure.main_body__7.timoschenko_input.set = [1, 2]
    # delete aerodynamic drag if requested
    if not withdrag:
        del htc.aerodrag
    # wind-speed  values
    htc.wind.tint = 0  # no TI
    htc.wind.turb_format = 0  # no turbulence
    htc.wind.tower_shadow_method = 0  # no tower shadow
    htc.wind.wsp = wsp  # mean wind speed
    htc.wind.shear_format = [1, wsp]  # constant wsp profile with height
    # define the append name based on the subfolder and the wind speed
    wsp_str = ('%.1f' % (wsp)).zfill(4)  # e.g., '05.0'
    if subfolder == '':
        append = f'_steady_{wsp_str}'  # e.g., '_05.0'
    else:
        append = f'_steady_{subfolder}_{wsp_str}'  # e.g., '_notiltnodragrigid_05.0'
    # update name and save file (reprint of _update_name_and_save() b.c. missing kwargs to set_name)
    save_dir = Path(htc_dir)  # sanitize inputs
    # get new name (excl extension) from HTCFile attribute "filename"
    name = Path(htc.filename).name.replace('.htc', append)
    # set filename using HTCFile method
    htc.set_name(name, resdir=res_dir, subfolder=subfolder, htcdir=htc_dir)
    # save the file
    htc.save((save_dir / subfolder / (name + '.htc')).as_posix())
    return


def main():
    """Create the htc files for the different cases, adjusting settings.
    Save the htc files in subfolders corresponding to the different cases.
    This code would be better placed at the end of your make_htc_files.py script...
    """
    # constants for this script
    del_htc_dir = True  # delete htc directory if it already exists?
    #master_htc = './_master/group7_3B_design.htc'
    master_htc = './hawc_files/our_design/htc/group7_3B_design_A3_part3_C1.htc' #this is maybe not the correct .htc file? - should it be one including contrller tuning?
    opt_path = './hawc_files/our_design/data/group7_3B_design_flex.opt'
    #cases = ['tilt', 'notilt', 'notiltrigid', 'notiltnodragrigid']
    cases = ['tilt', 'notilt']
    wsps = range(5, 25)  # wind speed range
    htc_dir = './htc_steady/'  # top-level folder to save htc files (can be path to gbar!)
    res_dir = './res_steady/'  # where HAWC2 should save res files, relative to its working directory
    # delete the top-level directory if requested
    _clean_directory(htc_dir, del_htc_dir)
    # make the files
    for case in cases:
        tilt = None  # default: don't change tilt
        rigid = False  # default: flexible blades and tower
        withdrag = True
        if 'notilt' in case:
            tilt = 0
        if 'rigid' in case:
            rigid = True
        if 'nodrag' in case:
            withdrag = False
        # generate the files
        for wsp in wsps:
            htc = MyHTC(master_htc)
            make_single_steady(htc, wsp,
                               htc_dir=htc_dir, res_dir=res_dir,
                               subfolder=case,
                               opt_path=opt_path,
                               tilt=tilt, rigid=rigid, withdrag=withdrag)


# the "script" part of this file
if __name__ == '__main__':
    main()
