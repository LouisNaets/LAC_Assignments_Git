"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/individual_design/_master/individual_design.htc'
    SAVE_HAWC2S_DIR = './hawc_files/individual_design'

    # make rigid hawc2s file for single-wsp opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_hawc2s_1wsp',
                    opt_path='./data/individual_design_1wsp.opt',
                    compute_steady_states=True,
                    save_power=True,
                    save_induction=True,
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))

    # make rigid hawc2s file for multi-tsr opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_hawc2s_multitsr',
                    opt_path='./data/individual_design_multitsr.opt',
                    compute_steady_states=True,
                    save_power=True,
                    save_induction=True,
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))
    
    # make rigid hawc2s file for multi-tsr 2nd opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_hawc2s_multitsr_new',
                    opt_path='./data/individual_design_multitsr_new.opt',
                    compute_steady_states=True,
                    save_power=True,
                    save_induction=True,
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))
    
    # make rigid hawc2s file for multi-tsr 3rd opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_hawc2s_multitsr_precise',
                    opt_path='./data/individual_design_multitsr_precise.opt',
                    compute_steady_states=True,
                    save_power=True,
                    save_induction=True,
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))
    
    # make rigid hawc2s file for new opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_compute_rigid_opt',
                    opt_path='./data/dtu_10mw_rigid.opt',
                    compute_optimal_pitch_angle=True,
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))
    
    # make flexible hawc2s file for new opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_compute_flex_opt',
                    opt_path='./data/dtu_10mw_rigid.opt',
                    compute_optimal_pitch_angle=True,
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))
    
    # make flexible hawc2s file (for A2)
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_flex',
                    opt_path='./data/individual_design_flex.opt',
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(0, 411.22))
    
    # make flexible (including min. rotational speed) hawc2s file (for A2)
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_flex_minrotspd',
                    opt_path='./data/individual_design_flex_minrotspd.opt',
                    gradient=False,
                    minpitch=0,
                    opt_lambda=7.15,
                    genspeed=(300, 411.22))