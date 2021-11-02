import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
import numpy as np
import os
import json

from grid_shape_lib.utils import diff_spec as diff_spec
from grid_shape_lib.utils import binning as binning

from grid_shape_lib.modules import grids as grids
from grid_shape_lib.utils import utils_analysis as ua

font = {'family': 'normal',
        'size': 12}
matplotlib.rc('font', **font)

# matplotlib.rc('xtick', labelsize=15)
# matplotlib.rc('ytick', labelsize=15)

SYM_LIST = 5 * ['.', 'o', 'v', '*', 's']
MYC = ['0', '0.20', '0.4', '0.6', '0.8']


class Layout:
    def __init__(
        self,
        path,
        pos,
        mask,
        mask_name,
        threshold,
        n_trig_thres,
        input_n_ring,
        primary="Proton"
    ):

        # "Creating" the layout
        # pos, _ = grids.create_grid_univ("trihex", 125, do_prune = False, input_n_ring = 10)
        self.pos = pos
        self.mask = mask
        self.mask_name = mask_name
        self.path = path
        self.threshold = threshold
        self.n_trig_thres = n_trig_thres
        self.plot_path = os.path.join(path, "plots")  # This is the directory where the plots will be saved
        self.events_data_dir = os.path.join(path, "events")  # This is the directory where "event" will be saved
        self.sanity_plots_dir = os.path.join(self.plot_path, "sanity_plots")
        self.primary = primary

        self.plot_suffix = "{}_{}_{}_{}".format(self.mask_name, self.primary, self.threshold, self.n_trig_thres)
        os.makedirs(self.plot_path, exist_ok=True)
        os.makedirs(self.events_data_dir, exist_ok=True)
        os.makedirs(self.sanity_plots_dir, exist_ok=True)

        self.merged_file_dir = self.path
        self.config_json_file = os.path.join(self.merged_file_dir, 'merge_config.json')

        self.xmin = pos[0].min()
        self.xmax = pos[0].max()

        self.ymin = pos[1].min()
        self.ymax = pos[1].max()

        with open(self.config_json_file) as f:
            self.config_merged = json.load(f)    ## read the merge config file. In this example it will not be used

        trihex_steps = self.config_merged["layouts"]["trihex"] ## here it is only 125, but can return a list of steps
        
        grid_shape = "trihex"
        self.input_n_ring = input_n_ring
        self.nb_hex = grids.hx.get_nb_hexagons(self.input_n_ring)
        self.area = np.int32(trihex_steps[0])**2 *3*np.sqrt(3)/2 * self.nb_hex  ## only valid for n_ring = 10

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
                primary,
                step, 
                self.threshold,
                self.n_trig_thres,
                prune_layout=(self.mask_name, self.mask),
                input_n_ring=input_n_ring
            )
            for step in trihex_steps
        ]

        ## This load the previously created arrays
        ev_select = [
            ua.get_ev_select(
                self.events_data_dir,
                "trihex",
                primary,
                step,
                self.threshold,
                self.n_trig_thres,
                prune_layout=(self.mask_name, self.mask)
            )
            for step in trihex_steps
        ]
        self.ev_select = np.concatenate([*ev_select])

        self.compute_energy_bins()
        self.compute_zenith_bins()

        self.delta_theta = self.zenith_bins_limits[1:] - \
            self.zenith_bins_limits[:-1]

        # self.delta_omega = 2*np.pi * self.delta_theta *np.pi/180 *\
        #    np.sin(np.pi/2 - self.zenith_bins_centers*np.pi/180)

        self.delta_omega = 2 * np.pi * (
            np.cos(np.pi/180*self.zenith_bins_limits[:-1]) - np.cos(np.pi/180*self.zenith_bins_limits[1:])
        ) 

        self.compute_trig_efficiency()
        self.compute_detection_rate()

    def compute_energy_bins(
        self,
        low_energy=16.3,
        high_energy=18.7,
        delta_log10=0.1
    ):
        ebl, ebc, de = binning.compute_energy_bins(low_energy, high_energy, delta_log10)

        self.energy_bins_limits = ebl
        self.energy_bins_centers = ebc
        self.delta_energy = de

    def compute_zenith_bins(self):

        zbl, zbc = binning.compute_zenith_bins(
            sec_theta_min=1.162,
            sec_theta_max=21.54,
            n_zenith_bins=16
        )

        self.zenith_bins_limits = zbl
        self.zenith_bins_centers = zbc

    def plot_layout(self, fig=1):
        # plot on the fly the pruned layout (i.e. plot is not saved)
        plt.figure(fig, figsize=(8, 6))
        plt.clf()
        plt.scatter(
            self.pos[0],
            self.pos[1],
            c=1-self.mask[:, 0],
            s=3
        )
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.axis('equal')
        plt.savefig(
            os.path.join(self.plot_path, "layout_{}.png".format(self.plot_suffix))
        )

    def plot_2d_efficiency(self, n_ener, n_zen, color='C0'):
        n_energy_bins = len(self.energy_bins_limits) - 1
        n_zenith_bins = len(self.zenith_bins_limits) - 1

        fig, axs = plt.subplots(
            n_zen,
            n_ener, 
            figsize=(10, 10),
            constrained_layout=True, 
            sharex='all', sharey='all',
    
        )
        fig1, axs1 = plt.subplots(
            n_zen,
            n_ener,   
            figsize=(10, 10),
            constrained_layout=True, 
            sharex='all', sharey='all',
        )
        fig2, axs2 = plt.subplots(
            n_zen,
            n_ener, 
            figsize=(10, 10),
            constrained_layout=True, 
            sharex='all', sharey='all',
        )
        fig3, axs3 = plt.subplots(
            n_zen,
            n_ener, 
            figsize=(10, 10),
            constrained_layout=True, 
            sharex='all', sharey='all',
        )

        ratio_ener = n_energy_bins//n_ener
        ratio_zen = n_zenith_bins//n_zen
   
        for i_ener in range(n_ener):
            for i_zen in range(n_zen):
                ind = np.where(
                    (self.ev_select[:,1] > self.energy_bins_limits[i_ener*ratio_ener]) * 
                    (self.ev_select[:,1] <= self.energy_bins_limits[(i_ener+1)*ratio_ener]) *
                    (self.ev_select[:,3] > self.zenith_bins_limits[i_zen*ratio_zen]) *
                    (self.ev_select[:,3] <= self.zenith_bins_limits[(i_zen+1)*ratio_zen]) 
                )
                ind2 = np.where( 
                    (self.ev_select[:,1] > self.energy_bins_limits[i_ener*ratio_ener]) * 
                    (self.ev_select[:,1] <= self.energy_bins_limits[(i_ener+1)*ratio_ener]) *
                    (self.ev_select[:,3] > self.zenith_bins_limits[i_zen*ratio_zen]) *
                    (self.ev_select[:,3] <= self.zenith_bins_limits[(i_zen+1)*ratio_zen]) *
                    (self.ev_select[:,4]==1)
                )
        
                vmax_denom = self.ev_select.shape[0] / n_ener / n_zen / (10-2) /(10-2)
                plt.figure(2008)

                x = self.ev_select[ind2[0], 6]
                y = self.ev_select[ind2[0], 7]
                
                x1 = self.ev_select[ind[0], 6]
                y1 = self.ev_select[ind[0], 7]
                
                theta = self.ev_select[ind2[0], 8] / 180 * np.pi
                theta1 = self.ev_select[ind[0], 8] / 180 * np.pi

                xp = x * np.cos(theta) - y * np.sin(theta)
                yp = x * np.sin(theta) + y * np.cos(theta)

                xp1 = x1 * np.cos(theta1) - y1 * np.sin(theta1)
                yp1 = x1 * np.sin(theta1) + y1 * np.cos(theta1)

                h_num = plt.hist2d(
                    x, y, bins=10, range = [[self.xmin, self.xmax],[self.ymin, self.ymax]]             
                )

                h_denom = plt.hist2d(
                    x1, y1, bins=10, range = [[self.xmin, self.xmax],[self.ymin, self.ymax]]
                )

                h_dist_triggered = plt.hist(np.sqrt(x**2 + y**2), bins=30, range = [0, np.sqrt(self.xmax**2 + self.ymax**2)])
                h_dist_all = plt.hist(np.sqrt(x1**2 + y1**2), bins=30, range = [0, np.sqrt(self.xmax**2 + self.ymax**2)])


                ax = axs[n_zen - 1 - i_zen, i_ener]
                ax1 = axs1[n_zen - 1 - i_zen, i_ener]
                ax2 = axs2[n_zen - 1 - i_zen, i_ener]
                ax3 = axs3[n_zen - 1 - i_zen, i_ener]

                pc = ax.pcolor(h_num[1], h_num[2], (h_num[0]/h_denom[0]).transpose(), vmin=0, vmax=1)
                pc1 = ax1.pcolor(h_num[1], h_num[2], h_num[0], vmin=0, vmax=vmax_denom)
                pc2 = ax2.pcolor(h_num[1], h_num[2], h_denom[0], vmin=0, vmax=vmax_denom)

                x_ = 0.5 * (h_dist_all[1][1:] + h_dist_all[1][:-1])
                pc3 = ax3.plot(x_, h_dist_triggered[0]/h_dist_all[0], color)
                ax3.set_ylim([0,1])
                ax.tick_params(labelsize=8)
                ax1.tick_params(labelsize=8)
                ax2.tick_params(labelsize=8)
                ax3.tick_params(labelsize=8)
                #ax.set_xticks(h_num[1])
                #ax.set_yticks(h_num[2])


#                axs1[i_zen, i_ener].pcolor(h_num[0], vmin=0, vmax=20)
#                axs2[i_zen, i_ener].pcolor(h_denom[0], vmin=0, vmax=40)
                if i_ener != 0:
                    ax.tick_params(labelleft=False)
                    ax1.tick_params(labelleft=False)
                    ax2.tick_params(labelleft=False)
                    ax3.tick_params(labelleft=False)
                if i_zen == n_zen - 1:
                    ax.set_title(
                        '%4.2f>log E/eV>%4.2f'%(
                            np.log10(1e18*self.energy_bins_limits[i_ener*ratio_ener]),
                            np.log10(1e18*self.energy_bins_limits[(i_ener+1)*ratio_ener])
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                    ax1.set_title(
                        '%4.2f>log E/eV>%4.2f'%(
                            np.log10(1e18*self.energy_bins_limits[i_ener*ratio_ener]),
                            np.log10(1e18*self.energy_bins_limits[(i_ener+1)*ratio_ener])
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                    ax2.set_title(
                        '%4.2f>log E/eV>%4.2f'%(
                            np.log10(1e18*self.energy_bins_limits[i_ener*ratio_ener]),
                            np.log10(1e18*self.energy_bins_limits[(i_ener+1)*ratio_ener])
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                    ax3.set_title(
                        '%4.2f>log E/eV>%4.2f'%(
                            np.log10(1e18*self.energy_bins_limits[i_ener*ratio_ener]),
                            np.log10(1e18*self.energy_bins_limits[(i_ener+1)*ratio_ener])
                        ),
                        fontdict={'fontsize': 'small'}
                    )



                if i_ener == n_ener - 1:
                    axs_ = ax.twinx()
                    axs1_ = ax1.twinx()
                    axs2_ = ax2.twinx()
                    axs3_ = ax3.twinx()
                    axs_.set_ylabel(
                        '%4.0f>zen>%4.0fdeg'%(
                            self.zenith_bins_limits[i_zen*ratio_zen],
                            self.zenith_bins_limits[(i_zen+1)*ratio_zen]
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                    axs1_.set_ylabel(
                        '%4.0f>zen>%4.0fdeg'%(
                            self.zenith_bins_limits[i_zen*ratio_zen],
                            self.zenith_bins_limits[(i_zen+1)*ratio_zen]
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                    axs2_.set_ylabel(
                        '%4.0f>zen>%4.0fdeg'%(
                            self.zenith_bins_limits[i_zen*ratio_zen],
                            self.zenith_bins_limits[(i_zen+1)*ratio_zen]
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                    axs3_.set_ylabel(
                        '%4.0f>zen>%4.0fdeg'%(
                            self.zenith_bins_limits[i_zen*ratio_zen],
                            self.zenith_bins_limits[(i_zen+1)*ratio_zen]
                        ),
                        fontdict={'fontsize': 'small'}
                    )
                   
                    axs_.tick_params(right=False, labelright=False)
                    axs1_.tick_params(right=False, labelright=False)
                    axs2_.tick_params(right=False, labelright=False)
                    axs3_.tick_params(right=False, labelright=False)

        #fig.tight_layout()
        #plt.subplots_adjust(wspace=0., 
        #           hspace=0.)
    
        fig.colorbar(pc, ax = axs[:, n_ener-1])
        fig.savefig(
            os.path.join(self.plot_path, "heatmap_{}.png".format(self.plot_suffix)),
            transparent=True
        )
        fig1.colorbar(pc1, ax = axs1[:, n_ener-1])
        fig1.savefig(
            os.path.join(self.plot_path, "heatmap_numerator{}.png".format(self.plot_suffix)),
            transparent=True
        )
        fig2.colorbar(pc2, ax = axs2[:, n_ener-1])
        fig2.savefig(
            os.path.join(self.plot_path, "heatmap_denominator{}.png".format(self.plot_suffix)),
            transparent=True
        )
        fig3.savefig(
            os.path.join(self.plot_path, "efficiency_array_{}.png".format(self.plot_suffix)),
            transparent=True
        )


    def compute_trig_efficiency(self):
        """
        compute the mean of variance of number of triggered antennas
        in various energy and zenith bins

        """
    
        n_energy_bins = len(self.energy_bins_limits) - 1
        n_zenith_bins = len(self.zenith_bins_limits) - 1

        mean_n_trig = np.zeros((n_energy_bins, n_zenith_bins))
        var_n_trig = np.zeros((n_energy_bins, n_zenith_bins))
        trig_efficiency = np.zeros((n_energy_bins, n_zenith_bins))

        


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
                    trig_efficiency[i_ener, i_zen] = 0
                else:
                    trig_efficiency[i_ener, i_zen] = sum(
                        self.ev_select[ind[0],4]
                    ) / np.size(self.ev_select[ind[0],4])



        self.mean_n_trig = mean_n_trig
        self.var_n_trig = var_n_trig
        self.trig_efficiency = trig_efficiency


    def compute_detection_rate(self):
        """
        compute the detection rate
        number of events per day per m2
        """
        n_energy_bins = len(self.energy_bins_centers)
        n_zenith_bins = len(self.zenith_bins_centers)

        detection_rate_old = np.zeros((n_energy_bins, n_zenith_bins))
        detection_rate = np.zeros((n_energy_bins, n_zenith_bins))

        differential_rate = np.zeros((n_energy_bins, n_zenith_bins))


        for i_ener, ener in enumerate(self.energy_bins_centers):
            for i_zen, zen in enumerate(self.zenith_bins_centers):
                differential_rate[i_ener, i_zen] = (
                        self.trig_efficiency[i_ener,i_zen] * 
                        diff_spec.tale_diff_flux(ener*1e18) * 
                        self.area * np.cos(zen*np.pi/180) / self.area*1e6 *
                        3600 * 24 # !! finally we divide by the area here
                    )


                detection_rate_old[i_ener, i_zen] = (
                        self.trig_efficiency[i_ener,i_zen] * 
                        diff_spec.tale_diff_flux(ener*1e18) * 
                        self.delta_energy[i_ener] *1e18 * self.delta_omega[i_zen] *
                        self.area * np.cos(zen*np.pi/180) / self.area*1e6 *
                        3600 * 24 # !! finally we divide by the area here
                    )

                dz = (self.zenith_bins_limits[1:] - self.zenith_bins_limits[:-1]) /180*np.pi
                detection_rate[i_ener, i_zen] = (
                        self.trig_efficiency[i_ener,i_zen] * 
                        diff_spec.tale_diff_flux(ener*1e18) * 
                        self.delta_energy[i_ener] *1e18 * np.cos(zen/180*np.pi) * 2*np.pi * dz[i_zen] *
                        self.area * np.cos(zen*np.pi/180) / self.area*1e6 *
                        3600 * 24 # !! finally we divide by the area here
                    )

        self.detection_rate = detection_rate
        self.detection_rate_old = detection_rate_old
        self.differential_rate = differential_rate
        

    def make_diff_event_rate(self, area_detector, do_plot=False):
        """
        Make the plot of differential event rate in day^-1 PeV^-1
        area_detector in km^2
        """
        n_energy_bins = len(self.energy_bins_centers)
        n_zenith_bins = len(self.zenith_bins_centers)

        diff_ev_rate = np.zeros((n_energy_bins, n_zenith_bins))
        diff_ev_rate_notrig = np.zeros((n_energy_bins, n_zenith_bins))

        diff_ev_rate = np.zeros((n_energy_bins, n_zenith_bins))
        diff_ev_rate_notrig = np.zeros((n_energy_bins, n_zenith_bins))

        for i_ener, ener in enumerate(self.energy_bins_centers):
            for i_zen, zen in enumerate(self.zenith_bins_centers):

                dz = (self.zenith_bins_limits[1:] - self.zenith_bins_limits[:-1]) /180*np.pi
                
                diff_ev_rate[i_ener, i_zen] = (
                        self.trig_efficiency[i_ener,i_zen]**1 * 
                        diff_spec.tale_diff_flux(ener*1e18) * 
                        np.cos(zen/180*np.pi) * 2*np.pi * dz[i_zen] *
                        np.cos(zen*np.pi/180) *
                        3600 * 24 * area_detector * 1e6 # !! finally we divide by the area here
                    )
                diff_ev_rate_notrig[i_ener, i_zen] = (
                        self.trig_efficiency[i_ener,i_zen]**0 * 
                        diff_spec.tale_diff_flux(ener*1e18) * 
                        np.cos(zen/180*np.pi) * 2*np.pi * dz[i_zen] *
                        np.cos(zen*np.pi/180) *
                        3600 * 24 * area_detector * 1e6 # !! finally we divide by the area here
                    )

        integrated_ev_rate = np.sum(diff_ev_rate, axis=1) * 1e15
        integrated_ev_rate_no_trig = np.sum(diff_ev_rate_notrig, axis=1) * 1e15
        
        self.integrated_ev_rate = integrated_ev_rate
        self.integrated_ev_rate_no_trig = integrated_ev_rate_no_trig
        
        if do_plot:
            plt.figure()
            plt.plot(self.energy_bins_centers, self.integrated_ev_rate)
            plt.plot(self.energy_bins_centers, self.integrated_ev_rate_no_trig)
            plt.ylabel('number of events [day'+"$^{-1}$"+"PeV"+"$^{-1}$"+']')
            plt.xlabel('Enegy [EeV]')
            plt.xscale('log')
            plt.yscale('log')


    def plot_trig_efficiency(self):
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
        trig_efficiency = self.trig_efficiency

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
                trig_efficiency[:,izen],
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
                trig_efficiency[i_bin, :],
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


    def plot_2D_differential_rate(self):

        fig, ax = plt.subplots()
        plt.clf()
        CS = plt.pcolor(
            np.log10(1e18*self.energy_bins_limits),
            self.zenith_bins_limits,
            self.differential_rate.T * 1e15,
            norm=LogNorm(), vmin=1e-7, vmax=1e2  
        )
        x = np.log10(1e18*self.energy_bins_centers)
        y = self.zenith_bins_centers
        xx,yy = np.meshgrid(x,y)
        levels = [1e-5, 1e-4, 1e-3, 1e-2]
        cs2 = plt.contour(xx,yy,
            self.differential_rate.T * 1e15, levels=levels, colors = 'k'
        )
        ax.clabel(cs2, inline=1, fontsize=10,colors="k", fmt= '%1.1e' )
        plt.xlabel('log E/eV')
        plt.ylabel('Zenith [deg]')
        plt.title(
            'differential event rate over array '+'$\\nu_{ev}$'+', {}'.format(
                self.plot_suffix
            )
        )
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('$ [day^{-1} km^{-2} PeV^{-1} sr^{-1}]$')
        plt.savefig(
            os.path.join(self.plot_path, "diff_ev_rate_{}.png".format(self.plot_suffix))
        )




    def plot_2D_trig_efficiency(self):
        fig, ax = plt.subplots()
        plt.clf()
        CS = plt.pcolor(
            np.log10(1e18*self.energy_bins_limits),
            self.zenith_bins_limits,
            self.trig_efficiency.T,
            norm=LogNorm(), vmin=1e-2, vmax=1
        )
        x = np.log10(1e18*self.energy_bins_centers)
        y = self.zenith_bins_centers
        xx,yy = np.meshgrid(x,y)
        levels = [0.01, 0.1, 0.5, 0.99]
        cs2 = plt.contour(xx,yy,
            self.trig_efficiency.T, levels=levels, colors = 'k'
        )
        ax.clabel(cs2, inline=1, fontsize=10,colors="k", fmt= '%1.2f' )
        
        plt.xlabel('log E/eV')
        plt.ylabel('Zenith [deg]')
        plt.title(
            'trig. efficiency, {}'.format(
                self.plot_suffix
            )
        )

        cbar = plt.colorbar(CS)
        #cbar.ax.set_ylabel('$$')
        plt.savefig(
            os.path.join(self.plot_path, "trig_efficiency_rate_{}.png".format(self.plot_suffix))
        )



    def plot_2D_detection_rate(self):
        fig, ax = plt.subplots()
        plt.clf()
        CS = plt.pcolor(
            np.log10(1e18*self.energy_bins_limits),
            self.zenith_bins_limits,
            self.detection_rate.T,
            norm=LogNorm(), vmin=1e-7, vmax=1e2
        )
        x = np.log10(1e18*self.energy_bins_centers)
        y = self.zenith_bins_centers
        xx,yy = np.meshgrid(x,y)
        levels = [0.01, 0.1, 1, 10]
        cs2 = plt.contour(xx,yy,
            self.detection_rate.T, levels=levels, colors = 'k'
        )
        ax.clabel(cs2, inline=1, fontsize=10,colors="k", fmt= '%1.2f' )
        
        plt.xlabel('log E/eV')
        plt.ylabel('Zenith [deg]')
        plt.title(
            'triggered event rate over array '+'$\\nu_{ev}$'+', {}'.format(
                self.plot_suffix
            )
        )

        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('$ [day^{-1} km^{-2}]$')
        plt.savefig(
            os.path.join(self.plot_path, "trig_ev_rate_{}.png".format(self.plot_suffix))
        )


    def plot_mean_n_trig(self):
        fig, ax = plt.subplots()
        CS = plt.pcolor(
            np.log10(1e18*self.energy_bins_limits),
            self.zenith_bins_limits,
            self.mean_n_trig.T,
            norm=LogNorm(), vmin=1, vmax=np.nanmax(self.mean_n_trig) 
        )
        x = np.log10(1e18*self.energy_bins_centers)
        y = self.zenith_bins_centers
        xx,yy = np.meshgrid(x,y)
        levels = [5, 10, 100]
        cs2 = plt.contour(xx,yy,
            self.mean_n_trig.T, levels=levels, colors = 'k'
        )
        ax.clabel(cs2, inline=1, fontsize=10,colors="k", fmt= '%2d' )
        plt.xlabel('log E/eV')
        plt.ylabel('Zenith [deg]')
        plt.title('Mean N trig, {}'.format(self.plot_suffix))
        cbar = plt.colorbar(CS)
        #cbar.ax.set_ylabel('$ [day^{-1} km^{-2} PeV^{-1} sr^{-1}]$')
        plt.savefig(
            os.path.join(
                self.plot_path,
                "mean_n_trig_{}.png".format(self.plot_suffix)
            )
        )


    def make_all_plots(self):
        self.plot_2D_detection_rate()
        self.plot_2D_differential_rate()
        self.plot_mean_n_trig()
        self.plot_2D_trig_efficiency()
        self.plot_layout()
        
