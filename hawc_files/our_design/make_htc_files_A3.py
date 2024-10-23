"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC
from pathlib import Path
from lacbox.io import load_ctrl_txt

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design.htc'
    SAVE_HAWC2S_DIR = './hawc_files/our_design'

    # load ctrl tuning data to dictionary
    fname = './hawc_files/our_design/res_hawc2s/group7_3B_design_A3_part2_C1_ctrl_tuning.txt'
    ctrltune_dict = load_ctrl_txt(fname)

    print('DICTIONARY KEYS:\n---------------------')
    [print(s) for s in ctrltune_dict.keys()]
    print('\nAERO_GAINS KEYS:\n---------------------')
    [print(s) for s in ctrltune_dict['aero_gains'].columns]
    
    # print the keys
    """  print('DICTIONARY KEYS:\n---------------------')
    [print(s) for s in ctrltune_dict.keys()]
    print('\nAERO_GAINS KEYS:\n---------------------')
    [print(s) for s in ctrltune_dict['aero_gains'].columns];
    DICTIONARY KEYS:
    ---------------------
    K_Nm/(rad/s)^2
    Irotor_kg*m^2
    KpTrq_Nm/(rad/s)
    KiTrq_Nm/rad
    KpPit_rad/(rad/s)
    KiPit_rad/rad
    K1_deg
    K2_deg^2
    Kp2_rad/(rad/s)
    Ko1_deg
    Ko2_deg^2
    aero_gains

    AERO_GAINS KEYS:
    ---------------------
    dq/dtheta_kNm/deg
    fit_kNm/deg
    dq/domega_kNm/(rad/s)
    fit_kNm/(rad/s) """
    

    # make htc file for tuning controller parameters
    htc = MyHTC(ORIG_PATH)
    print(dir(htc))  # Check if make_hawc2s_ctrltune exists

    htc.make_hawc2s_ctrltune(SAVE_HAWC2S_DIR,
                    rigid=False,
                    gradient = True,
                    append='_controller_tuning',
                    opt_path='./data/group7_3B_design_flex.opt',
                    opt_lambda=7.5,
                    partial_load=(0.05, 0.7),
                    full_load=(0.06, 0.7),
                    compute_steady_states=True,
                    save_power=True,
                    compute_controller_input=True)

    omega_Omegas = [0.05, 0.01, 0.1]
    for idx, omega in enumerate(omega_Omegas, start=1):
        htc = MyHTC(ORIG_PATH)
        append_str = f'_A3_part2_C{idx}'
        htc.make_hawc2s_ctrltune(SAVE_HAWC2S_DIR,
                    rigid=False,
                    gradient = True,
                    append=append_str,
                    opt_path='./data/group7_3B_design_flex.opt',
                    opt_lambda=7.5,
                    constant_power=1,
                    full_load=(omega, 0.7),
                    compute_steady_states=True,
                    save_power=True,
                    compute_controller_input=True)
        
    for idx, omega in enumerate(omega_Omegas, start=4):
        htc = MyHTC(ORIG_PATH)
        append_str = f'_A3_part2_C{idx}'
        htc.make_hawc2s_ctrltune(SAVE_HAWC2S_DIR,
                    rigid=False,
                    gradient = True,
                    append=append_str,
                    opt_path='./data/group7_3B_design_flex.opt',
                    opt_lambda=7.5,
                    constant_power=0,
                    full_load=(omega, 0.7),
                    compute_steady_states=True,
                    save_power=True,
                    compute_controller_input=True)