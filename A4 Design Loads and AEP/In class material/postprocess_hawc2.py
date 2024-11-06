# -*- coding: utf-8 -*-
"""Post-process a selection of HAWC2 files run on the cluster.
"""
from pathlib import Path

from lacbox.postprocess import process_statistics


# inputs
res_dir_A = Path('./res_turb_2/tca')  # directory with res files to process
res_dir_B = Path('./res_turb_2/tcb')  # directory with res files to process
calc_del = False  # calculate DELs in the statistics? It takes longer.
save_path_A = './res_turb_2/tca/group7_turbA_stats.csv'  # where should I save the stats file?
save_path_B = './res_turb_2/tcb/group7_turbB_stats.csv'  # where should I save the stats file?

# call the function
stats_df = process_statistics(res_dir_A, save_path_A)
stats_df = process_statistics(res_dir_B, save_path_B)

