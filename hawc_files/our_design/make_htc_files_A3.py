"""Make all the htc files for the LAC course from a single base file.

Requires myteampack (which requires lacbox).
"""
from myteampack import MyHTC

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