"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC

if __name__ == '__main__':
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design.htc'
    SAVE_HAWC2S_DIR = './hawc_files/our_design'

    # make htc file for tuning controller parameters
    htc = MyHTC(ORIG_PATH)
    htc.make_hawc2s_ctrltune(save_dir=SAVE_HAWC2S_DIR,
                    rigid=False,
                    append='_controller_tuning',
                    opt_path='./data/dtu_10mw_rigid.opt')