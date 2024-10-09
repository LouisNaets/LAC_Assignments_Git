"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design_old_st.htc'
    SAVE_HAWC2S_DIR = './hawc_files/our_design'
    # make rigid hawc2s file for single-wsp opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_hawc2s_1wsp',
                    opt_path='./data/group7_3B_design_1wsp.opt',
                    compute_steady_states=True,
                    save_power=True,
                    save_induction=True)

    # make rigid hawc2s file for multi-tsr opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_hawc2s_multitsr',
                    opt_path='./data/group7_3B_design_multitsr.opt',
                    compute_steady_states=True,
                    save_power=True,
                    save_induction=True)
    
    # make rigid hawc2s file for new opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=True,
                    append='_compute_rigid_opt',
                    opt_path='./data/dtu_10mw_rigid.opt',
                    compute_optimal_pitch_angle=True,
                    minpitch=0,
                    opt_lambda=7.5,
                    genspeed=(0, 431.35))
    
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design.htc'
    # make flexible hawc2s file for new opt file
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_compute_flex_opt',
                    opt_path='./data/dtu_10mw_rigid.opt',
                    compute_optimal_pitch_angle=True,
                    minpitch=0,
                    opt_lambda=7.5,
                    genspeed=(0, 431.35))
    
    # make flexible hawc2s file (for A2)
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_flex',
                    opt_path='./data/group7_3B_design_flex.opt',
                    minpitch=0,
                    opt_lambda=7.5,
                    genspeed=(0, 431.35))