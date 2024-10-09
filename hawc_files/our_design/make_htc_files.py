"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design_old_st.htc'
    SAVE_HAWC2S_DIR = './hawc_files/our_design'

    # make htc file for tuning controller parameters
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s(SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_controller_tuning',
                    opt_path='./data/dtu_10mw_rigid.opt',
                    torque_controller={'natural_frequency': 0.05, 'damping_ratio': 0.7},
                    pitch_controller={'natural_frequency': 0.06, 'damping_ratio': 0.7})