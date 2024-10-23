"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC
from lacbox.io import load_ctrl_txt

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design.htc'
    SAVE_HAWC2S_DIR = './hawc_files/our_design'

    # make htc file for tuning controller parameters
    htc = MyHTC(ORIG_PATH)
    # print(dir(htc))  # Check if make_hawc2s_ctrltune exists

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
    
    #Part 3 requires a new subfolder for saving htc files
    SAVE_HAWC2S_DIR = './hawc_files/our_design/htc'
    cp_dict = load_ctrl_txt('./hawc_files/our_design/res_hawc2s/group7_3B_design_controller_tuning_ctrl_tuning.txt')
    htc = MyHTC(ORIG_PATH)
    htc.make_step(save_dir=SAVE_HAWC2S_DIR,
                  append="_A3_part3",
                  cp_dict=cp_dict, 
                  t_start=0.,
                  t_end=880.,
                  start_wsp=10.,
                  tint=0.,
                  turb_format=0, 
                  shear_format=(3,0),
                  tower_shadow_method=0,
                  wind_ramp_abs=(0, 880, 4, 25))
