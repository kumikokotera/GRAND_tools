import matplotlib.pyplot as plt
import numpy as np
import os
import json
import diff_spec
import grids as grids
import utils_analysis as ua





#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/InterpolationOutputExample/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/StshpLibrary02/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/P2PdataNew/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/StshpLibrary-HDF5-Grids/"
path = "/Users/benoitl/Documents/GRAND/StshpLibaryHDF5-Grids/"

#plot_path = '/Users/kotera/BROQUE/Plots_GRAND/'
plot_path = '/Users/benoitl/Documents/GRAND/plots_tests'
#events_data_dir = "/Users/kotera/BROQUE/Data_GRAND/Matias/event_data"
events_data_dir = "/Users/benoitl/Documents/GRAND/event_data"

os.makedirs(plot_path, exist_ok=True)
os.makedirs(events_data_dir, exist_ok=True)

threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 

merged_file_dir = "/Users/benoitl/Documents/GRAND/StshpLibrary-HDF5-Grids_merge/"

config_json_file = os.path.join(merged_file_dir, 'merge_config.json')

with open(config_json_file) as f:
    config_merged = json.load(f)




grid_shape = "trihex"
radius = 1000



pos, offset, mask = grids.create_grid_univ(grid_shape, radius, do_trim=True)


trihex_steps = config_merged["layouts"]["trihex"]
hexhex_steps = config_merged["layouts"]["hexhex"]


[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        grid_shape,
        "Proton",
        step, 
        threshold,
        n_trig_thres,
        prune_layout=("prune1", mask)
    )
    for step in trihex_steps
]



ev_select_hexhex = [
    ua.get_ev_select(
        events_data_dir,
        "hexhex",
        "Proton",
        step,
        threshold,
        n_trig_thres
    )
    for step in hexhex_steps
]
ev_select_hexhex = np.concatenate([*ev_select_hexhex])  


ev_select_trihex = [
    ua.get_ev_select(
        events_data_dir,
        "trihex",
        "Proton",
        step,
        threshold,
        n_trig_thres,
        prune_layout=("prune1", mask)
    )
    for step in trihex_steps
]
ev_select_trihex = np.concatenate([*ev_select_trihex])  











plt.figure(1)
plt.clf()
plt.scatter(pos[0,:], pos[1,:], c = mask[:,0])
plt.axis("equal")
plt.colorbar()






plt.show()