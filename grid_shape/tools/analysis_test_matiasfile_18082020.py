import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 


path = "/Users/benoitl/Documents/GRAND/test_matias_18082020/"

plot_path = os.path.join(path, "plots")
events_data_dir = os.path.join(path, "events")

os.makedirs(plot_path, exist_ok=True)
os.makedirs(events_data_dir, exist_ok=True)

threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 

merged_file_dir = path
config_json_file = os.path.join(merged_file_dir, 'merge_config.json')

with open(config_json_file) as f:
    config_merged = json.load(f)







#############################################################
## focusing in trihex nring = 10, step = 125
#############################################################



pos, offset = grids.create_grid_univ("trihex", 125, do_prune=False, input_n_ring=10)
# Make pruning mask
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 125, do_prune=True, input_n_ring=10)


plt.figure()
plt.scatter(pos2[0], pos2[1], c = mask2[:,0], s=1)


trihex_steps = config_merged["layouts"]["trihex"]

grid_shape = "trihex"

## Create ev_select with pruning 
[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        grid_shape,
        "Gamma",
        step, 
        threshold,
        n_trig_thres,
        prune_layout=("simple", mask2)
    )
    for step in trihex_steps
]



## Create ev_select without pruning 
[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        grid_shape,
        "Gamma",
        step, 
        threshold,
        n_trig_thres,
    )
    for step in trihex_steps
]



ev_select_trihex_pruning = [
    ua.get_ev_select(
        events_data_dir,
        "trihex",
        "Gamma",
        step,
        threshold,
        n_trig_thres,
        prune_layout=("simple", mask2)
    )
    for step in trihex_steps
]
ev_select_trihex_pruning= np.concatenate([*ev_select_trihex_pruning])  



ev_select_trihex_nopruning = [
    ua.get_ev_select(
        events_data_dir,
        "trihex",
        "Gamma",
        step,
        threshold,
        n_trig_thres
    )
    for step in trihex_steps
]
ev_select_trihex_nopruning= np.concatenate([*ev_select_trihex_nopruning])  



plt.figure()
plt.hist(ev_select_trihex_nopruning[:,0], bins = 500, range=[0, 1100], label = 'No pruning')
plt.hist(ev_select_trihex_pruning[:,0], bins = 500, range=[0, 1100], alpha=0.5, label = 'pruning')
plt.ylabel("# of event")
plt.xlabel('Num_triggered')
plt.legend(loc=0)







