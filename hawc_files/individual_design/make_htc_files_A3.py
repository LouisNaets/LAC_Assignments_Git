"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC
from lacbox.io import load_ctrl_txt
from sys import exit

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/individual_design/_master/individual_design.htc'
    SAVE_HAWC2S_DIR = './hawc_files/individual_design'

    # make htc file for tuning controller parameters
    htc = MyHTC(ORIG_PATH)
    print(dir(htc))  # Check if make_hawc2s_ctrltune exists

    htc.make_hawc2s_ctrltune(SAVE_HAWC2S_DIR,
                    rigid=False,
                    gradient = True,
                    append='_controller_tuning',
                    opt_path='./data/individual_design_flex_minrotspd.opt',
                    opt_lambda=7.15,
                    genspeed=(300, 411.22), #added these after the part 2 simulation
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
                    opt_path='./data/individual_design_flex_minrotspd.opt',
                    opt_lambda=7.15,
                    genspeed=(300, 411.22), #added these after the part 2 simulation
                    constant_power=1,
                    full_load=(omega, 0.7),
                    compute_steady_states=True,
                    save_power=True,
                    compute_controller_input=True)
        
    test = False
    if test:
        for idx, omega in enumerate(omega_Omegas, start=4):
            htc = MyHTC(ORIG_PATH)
            append_str = f'_A3_part2_C{idx}'
            htc.make_hawc2s_ctrltune(SAVE_HAWC2S_DIR,
                        rigid=False,
                        gradient = True,
                        append=append_str,
                        opt_path='./data/individual_design_flex_minrotspd.opt',
                        opt_lambda=7.15,
                        genspeed=(300, 411.22), #added these after the part 2 simulation
                        constant_power=0,
                        full_load=(omega, 0.7),
                        compute_steady_states=True,
                        save_power=True,
                        compute_controller_input=True)
    
    

    #Part 3 requires a new subfolder for saving htc files
    SAVE_HAWC2S_DIR = './hawc_files/individual_design/htc'
    cp_dict = load_ctrl_txt('./hawc_files/individual_design/res_hawc2s/individual_design_controller_tuning_ctrl_tuning.txt')
    htc = MyHTC(ORIG_PATH)
    htc.make_step(save_dir=SAVE_HAWC2S_DIR,
                  append="_A3_part3",
                  cp_dict=cp_dict, 
                  t_start=0.,
                  t_end=1862.,
                  start_wsp=10.,
                  tint=0.,
                  turb_format=0, 
                  shear_format=(3,0),
                  tower_shadow_method=0,
                  wind_ramp_abs=(0, 1862, 4, 25))
    
    test = False
    if test:
        for idx in range(1,7):
            htc = MyHTC(ORIG_PATH)
            fname = f'./hawc_files/individual_design/res_hawc2s/individual_design_A3_part2_C{idx}_ctrl_tuning.txt'
            ctrltune_dict = load_ctrl_txt(fname)
            #print('DICTIONARY KEYS:\n---------------------')
            #[print(s) for s in ctrltune_dict.keys()]
            append_str = f'_A3_part3_C{idx}'
            htc.make_step(save_dir=SAVE_HAWC2S_DIR,
                    append=append_str,
                    cp_dict=ctrltune_dict, 
                    t_start=0.,
                    t_end=1862., #this includes the 100s transient
                    start_wsp=4.,
                    tint=0.,
                    turb_format=0, 
                    shear_format=(3,0),
                    tower_shadow_method=0,
                    wind_ramp_abs=(0, 1862, 4, 25))
        
    C7_list = ['0.03_0.7', '0.03_0.8', '0.03_0.9', '0.04_0.8', '0.04_0.9', '0.05_0.8', '0.05_0.9']

    idx = 0
    for conditions in C7_list:
        idx = idx+1
        htc = MyHTC(ORIG_PATH)
        fname = f'./hawc_files/individual_design/res_hawc2s/individual_design_CNEW_{conditions}.txt'
        ctrltune_dict = load_ctrl_txt(fname)
        append_str = f'_CNEW_{idx}'
        htc.make_step(save_dir=SAVE_HAWC2S_DIR,
                  append=append_str,
                  cp_dict=ctrltune_dict, 
                  t_start=0.,
                  t_end=1862., #this includes the 100s transient
                  start_wsp=4.,
                  tint=0.,
                  turb_format=0, 
                  shear_format=(3,0),
                  tower_shadow_method=0,
                  wind_ramp_abs=(0, 1862, 4, 25))