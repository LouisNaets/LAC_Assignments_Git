"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC


if __name__ == '__main__':
    ORIG_PATH = './hawc_files/our_design/_master/group7_3B_design.htc'
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