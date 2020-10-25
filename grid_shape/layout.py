import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 
font = {'family' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

#matplotlib.rc('xtick', labelsize=15) 
#matplotlib.rc('ytick', labelsize=15) 



SYM_LIST = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
MYC = ['0','0.20','0.4','0.6','0.8']


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


        trihex_steps  =  self.config_merged["layouts"]["trihex"] ## here it is only 125, but can return a list of steps
        
        grid_shape  =  "trihex"

        self.area = np.int32(trihex_steps[0])**2 *3*np.sqrt(3)/2 * 331  ## only valid for n_ring = 10

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
   
        self.delta_theta = self.zenith_bins_limits[1:] - \
            self.zenith_bins_limits[:-1]

        self.delta_omega = 2*np.pi * self.delta_theta *np.pi/180 * np.sin(
            np.pi/2 - self.zenith_bins_centers*np.pi/180)

        self.compute_trig_rate()
        self.compute_detection_rate()
         

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
        n_zenith_bins = 16
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
        plt.figure(fig, figsize=(8, 6))
        plt.clf()
        plt.scatter(
            self.pos[0],
            self.pos[1],
            c = 1-self.mask[:,0],
            s = 3
        )
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.axis('equal')



    def compute_trig_rate(self):
        """
        compute the mean of variance of number of triggered antennas
        in various energy and zenith bins

        """
    
        n_energy_bins = len(self.energy_bins_limits) - 1
        n_zenith_bins = len(self.zenith_bins_limits) - 1

        mean_n_trig = np.zeros((n_energy_bins, n_zenith_bins))
        var_n_trig = np.zeros((n_energy_bins, n_zenith_bins))
        trig_rate = np.zeros((n_energy_bins, n_zenith_bins))


        for i_ener in range(n_energy_bins):
            for i_zen in range(n_zenith_bins):
                ind = np.where(
                    (self.ev_select[:,1] > self.energy_bins_limits[i_ener]) * 
                    (self.ev_select[:,1] <= self.energy_bins_limits[i_ener+1]) *
                    (self.ev_select[:,3] > self.zenith_bins_limits[i_zen]) *
                    (self.ev_select[:,3] <= self.zenith_bins_limits[i_zen+1])
                )

                mean_n_trig[i_ener, i_zen] = np.mean(self.ev_select[ind[0],0])
                var_n_trig[i_ener, i_zen] = np.var(self.ev_select[ind[0],0])
                if len(ind[0]) == 0:
                    trig_rate[i_ener, i_zen] = 0
                else:
                    trig_rate[i_ener, i_zen] = sum(
                        self.ev_select[ind[0],4]
                    ) / np.size(self.ev_select[ind[0],4])

        self.mean_n_trig = mean_n_trig
        self.var_n_trig = var_n_trig
        self.trig_rate = trig_rate


    def compute_detection_rate(self):
        """
        compute the detection rate
        number of events per day per m2
        """
        n_energy_bins = len(self.energy_bins_centers)
        n_zenith_bins = len(self.zenith_bins_centers)

        detection_rate = np.zeros((n_energy_bins, n_zenith_bins))
        
        for i_ener, ener in enumerate(self.energy_bins_centers):
            for i_zen, zen in enumerate(self.zenith_bins_centers):

                detection_rate[i_ener, i_zen] = (
                        self.trig_rate[i_ener,i_zen] * 
                        diff_spec.tale_diff_flux(ener*1e18) * 
                        self.delta_energy[i_ener] *1e18 * self.delta_omega[i_zen] *
                        self.area * np.cos(zen*np.pi/180) / self.area*1e6 # !! finally we divide by the area here
                    )
        self.detection_rate = detection_rate


    def plot_trig_rate(self):
        """
        plot Ntriggered antennas vs energies for fixed steps
        enerbins and zenbins are the arrays of the bin limits
        and plot the rate of triggering events
        """
        meanNtrig_ener = self.mean_n_trig
        varNtrig_ener = self.var_n_trig
        enerbins = self.energy_bins_limits
        enerbins_center = self.energy_bins_centers
        zenbins = self.zenith_bins_limits
        zenbins_center = self.zenith_bins_centers
        layout = self.mask_name
        trig_rate = self.trig_rate

        plt.figure()
        for izen in range(0, len(zenbins)-1):
            plt.errorbar(
                enerbins_center,
                meanNtrig_ener[:,izen],
                yerr=np.sqrt(varNtrig_ener[:,izen]), 
                fmt=SYM_LIST[izen],
                capsize=2,
                alpha=0.7,
                ms=7,
                label='%4.0f > zen >%4.0f deg'%(zenbins[izen],zenbins[izen+1])
            )
            
            plt.yscale('log')
            plt.ylabel('N triggered antennas')
            plt.xlabel('energy [EeV]')
            plt.xscale('log')
            plt.title('%s'%(layout))
            plt.legend(loc=4)

        plt.figure()
        n_energy_bins = len(enerbins_center)
        n_bins_to_plot = 5
        i_bins = [(i+1) * n_energy_bins // (n_bins_to_plot+1) for i in range(n_bins_to_plot)]
        for k, i_bin in enumerate(i_bins):
            plt.errorbar(
                zenbins_center,
                meanNtrig_ener[i_bin, :],
                yerr=np.sqrt(varNtrig_ener[i_bin, :]), 
                fmt=SYM_LIST[k],
                capsize=2,
                alpha=0.7,
                ms=7,
                label='%4.2f > log E/eV >%4.2f'%(
                    np.log10(enerbins[i_bin]*1e18),
                    np.log10(enerbins[i_bin+1]*1e18)
                )
            )

        plt.yscale('log')
        plt.ylabel('N triggered antennas')
        plt.xlabel('Zenith [deg]')
        plt.title('%s'%(layout))
        plt.legend(loc=0, ncol=2)


        plt.figure()
        for izen in range(0, len(zenbins)-1):
            plt.errorbar(
                enerbins_center,
                trig_rate[:,izen],
                fmt=SYM_LIST[izen],
                ls='-',
                capsize=2,
                alpha=0.7,
                ms=7,
                label='%4.0f > zen >%4.0f deg'%(zenbins[izen],zenbins[izen+1])
            )
            
            plt.yscale('log')
            plt.ylabel('Triggered event rate')
            plt.xlabel('energy [EeV]')
            plt.xscale('log')
            plt.title('%s'%(layout))
            plt.legend(loc=4)
    
        plt.figure()
        n_energy_bins = len(enerbins_center)
        n_bins_to_plot = 5
        i_bins = [(i+1) * n_energy_bins // (n_bins_to_plot+1) for i in range(n_bins_to_plot)]
        for k, i_bin in enumerate(i_bins):
            plt.errorbar(
                zenbins_center,
                trig_rate[i_bin, :],
                fmt=SYM_LIST[k],
                ls='-',
                capsize=2,
                alpha=0.7,
                ms=7,
                label='%4.2f > log E/eV >%4.2f'%(
                    np.log10(enerbins[i_bin]*1e18),
                    np.log10(enerbins[i_bin+1]*1e18)
                )
            )

        plt.yscale('log')
        plt.ylabel('Triggered event rate')
        plt.xlabel('Zenith [deg]')
        plt.title('%s'%(layout))
        plt.legend(loc=0, ncol=2)





 
