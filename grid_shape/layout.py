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
        #pos, _  =  grids.create_grid_univ("trihex", 125, do_prune = False, input_n_ring = 10)
        self.pos  =  pos
        self.mask  =  mask
        self.mask_name  =  mask_name
        self.path  =  path
        self.threshold  =  threshold
        self.n_trig_thres  =  n_trig_thres
        self.plot_path  =  os.path.join(path, "plots")  # This is the directory where the plots will be saved
        self.events_data_dir  =  os.path.join(path, "events") # This is the directory where "event" will be saved
        self.sanity_plots_dir  =  os.path.join(self.plot_path, "sanity_plots")

        os.makedirs(self.plot_path, exist_ok = True)
        os.makedirs(self.events_data_dir, exist_ok = True)
        os.makedirs(self.sanity_plots_dir, exist_ok = True)

        self.merged_file_dir  =  self.path
        self.config_json_file  =  os.path.join(self.merged_file_dir, 'merge_config.json')


        with open(self.config_json_file) as f:
            self.config_merged  =  json.load(f)    ## read the merge config file. In this example it will not be used


        trihex_steps  =  self.config_merged["layouts"]["trihex"]  ## here it is only 125, but can return a list of steps

        grid_shape  =  "trihex"

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
                prune_layout = (self.mask_name, self.mask), 
                input_n_ring = 10
            )
            for step in trihex_steps
        ]

        ## This load the previously created arrays
        ev_select =  [
            ua.get_ev_select(
                self.events_data_dir,
                "trihex",
                "Gamma",
                step,
                self.threshold,
                self.n_trig_thres,
                prune_layout = (self.mask_name, self.mask)
            )
            for step in trihex_steps
        ]
        self.ev_select  =  np.concatenate([*ev_select])  

        self.compute_energy_bins()
        self.compute_zenith_bins()

    def compute_energy_bins(self,
        low_energy = 16,
        high_energy = 19,
        delta_log10 = 0.1
    ):

        dum_e = np.arange(low_energy, high_energy+delta_log10/10, delta_log10) 
        self.energy_bins_limits  =  10**(
            dum_e - delta_log10/2 
        ) / 1e18

        self.energy_bins_centers =  10**(dum_e[0:-1])/1e18
        self.delta_energy  =  self.energy_bins_limits[1:] - self.energy_bins_limits[:-1]

    def compute_zenith_bins(self):

        max_theta = np.rad2deg(np.arccos(1.0/21.54))
        min_theta = np.rad2deg(np.arccos(1.0/1.162))
        sec_theta_min = 1.0/np.cos(np.deg2rad(min_theta))
        sec_theta_max = 1.0/np.cos(np.deg2rad(max_theta))
        n_zenith_bins = 8
        print("nbins:"+str(n_zenith_bins))
        print("max theta:"+str(max_theta)+" -> " + str(sec_theta_max))
        print("min theta:"+str(min_theta)+" -> " + str(sec_theta_min))
        log_sec_theta_min = np.log10(sec_theta_min)
        log_sec_theta_max = np.log10(sec_theta_max)
        log_increment = (log_sec_theta_max-log_sec_theta_min)/n_zenith_bins
        print("On log spacing")
        print("max log_sec_theta:"+str(max_theta)+" -> " + str(log_sec_theta_max))
        print("min log_sec_theta:"+str(min_theta)+" -> " + str(log_sec_theta_min))
        print("log increment is sec theta:"+str(log_increment)+" a factor " + str(np.power(10,log_increment)))
        limits_log_secant_bins = np.logspace(
            log_sec_theta_min,
            log_sec_theta_max,
            n_zenith_bins+1,
            base=10
        )
        print("log_limits:" + str(limits_log_secant_bins))
        center_log_secant_bins = limits_log_secant_bins[0:n_zenith_bins]*np.power(10,log_increment/2)
        print("log centers:"+ str(center_log_secant_bins))
        self.zenith_bins_limits = np.rad2deg(np.arccos(1.0/limits_log_secant_bins))
        print("limits:" + str(self.zenith_bins_limits))
        self.zenith_bins_centers = np.rad2deg(np.arccos(1.0/center_log_secant_bins))
        print("centers:" + str(self.zenith_bins_centers))
       

    def plot_layout(self, fig = 1):
        ## plot on the fly the pruned layout (i.e. plot is not saved)
        plt.figure(fig)
        plt.clf()
        plt.scatter(
            self.pos[0],
            self.pos[1],
            c = self.mask[:,0],
            s = 1
        )




# class SubLayout(Layout):
#     def __init__(self, )





 
