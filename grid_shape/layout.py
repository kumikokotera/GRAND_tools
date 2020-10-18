import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 

class Layout:
    def __init__(self, path, pos, mask, mask_name, threshold, n_trig_thres):

        # "Creating" the layout
        #pos, _ = grids.create_grid_univ("trihex", 125, do_prune=False, input_n_ring=10)
        self.pos = pos
        self.mask = mask
        self.mask_name = mask_name
        self.path = path
        self.threshold = threshold
        self.n_trig_thres = n_trig_thres
        self.plot_path = os.path.join(path, "plots")  # This is the directory where the plots will be saved
        self.events_data_dir = os.path.join(path, "events") # This is the directory where "event" will be saved
        self.sanity_plots_dir = os.path.join(self.plot_path, "sanity_plots")

        os.makedirs(self.plot_path, exist_ok=True)
        os.makedirs(self.events_data_dir, exist_ok=True)
        os.makedirs(self.sanity_plots_dir, exist_ok=True)

        self.merged_file_dir = self.path
        self.config_json_file = os.path.join(self.merged_file_dir, 'merge_config.json')

        with open(self.config_json_file) as f:
            self.config_merged = json.load(f)    ## read the merge config file. In this example it will not be used


        trihex_steps = self.config_merged["layouts"]["trihex"]  ## here it is only 125, but can return a list of steps

        grid_shape = "trihex"

        ## Create ev_select with pruning.
        ## The ev_select is the selection of the events that match the threshold and n_trig_thres criteria
        ## so the next line create an numpy array with:
        # - the number of triggered antenna (yes if Ep2p_tot > threshold),
        # - the energy of the event,
        # - the step
        # - the zenith
        # - and if the event triggers (number of triggred antenna > n_trig_thres )
        [
            ua.create_ev_select(
                self.events_data_dir,
                self.merged_file_dir,
                self.sanity_plots_dir,
                grid_shape,
                "Gamma",
                step, 
                self.threshold,
                self.n_trig_thres,
                prune_layout=(self.mask_name, self.mask), 
                input_n_ring=10
            )
            for step in trihex_steps
        ]

        ## This load the previously created arrays
        ev_select= [
            ua.get_ev_select(
                self.events_data_dir,
                "trihex",
                "Gamma",
                step,
                self.threshold,
                self.n_trig_thres,
                prune_layout=(self.mask_name, self.mask)
            )
            for step in trihex_steps
        ]
        self.ev_select = np.concatenate([*ev_select])  






    def plot_layout(self, fig=1):
        ## plot on the fly the pruned layout (i.e. plot is not saved)
        plt.figure(fig)
        plt.clf()
        plt.scatter(
            self.pos[0],
            self.pos[1],
            c=self.mask[:,0],
            s=1
        )




# class SubLayout(Layout):
#     def __init__(self, )





 
