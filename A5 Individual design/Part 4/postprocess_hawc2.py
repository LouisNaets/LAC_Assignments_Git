# -*- coding: utf-8 -*-
"""Post-process a selection of HAWC2 files run on the cluster.
"""
from pathlib import Path

from lacbox.postprocess import process_statistics


# inputs
res_dir = Path('./A5 Individual design/Part 4/tcb')  # directory with res files to process
calc_del = False  # calculate DELs in the statistics? It takes longer.
save_path = './A5 Individual design/Part 4/individual_design_turb_tcb_stats.csv'  # where should I save the stats file?

# call the function
stats_df = process_statistics(res_dir, save_path)
