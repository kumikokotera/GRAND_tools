import numpy as np
import os
import json
import matplotlib.pyplot as plt
import ijson
from grid_shape import grids as grids

'''
Definition of class to read the simulation output files
and of event selection routines 
to perform simulation output analysis in the analysis.py routine
'''

SYM_LIST = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
MYC = ['0','0.20','0.4','0.6','0.8']




class Event_old:
    def __init__(self,f1,f2, step, name):
        p2px, p2py, p2pz, p2ptot = np.loadtxt(f1)
        self.p2px = p2px
        self.p2py = p2py
        self.p2pz = p2pz
        self.p2ptot = p2ptot
        # JobName,Primary,Energy,Zenith,Azimuth,XmaxDistance,SlantXmax,RandomCore[0],RandomCore[1],RandomAzimuth,HadronicModel
        A = open(f2).readlines()[0]
        A = A.strip().split()
        self.jobname = A[0]
        self.primary = A[1]
        self.energy = np.float32(A[2])
        self.zenith = np.float32(A[3])
        self.azimuth = np.float32(A[4])
        self.xmax_distance = np.float32(A[5])
        self.slant_xmax = np.float32(A[6])
        self.random_core0 = np.float32(A[7])
        self.random_core1 = np.float32(A[8])
        self.random_azimuth = np.float32(A[9])
        self.hadronic_model = A[10]
        self.step = np.float32(step)
        self.name = name
        self.init_layout()

    def init_layout(self):
        if "hexhex" in self.name:
            self.layout = "hexhex"
        elif "rect" in self.name:
            self.layout = "rect"
        elif "trihex" in self.name:
            self.layout = "trihex"
        else:
            self.layout = "unknown"

    def is_triggered1(self, threshold):
        return self.p2ptot > threshold


class Event:
    def __init__(self, showerinfo, showerdata, step, name, prune_layout=()):
       
        p2px, p2py, p2pz, p2ptot = showerdata
        if prune_layout == ():
            self.p2px = np.array(p2px)
            self.p2py = np.array(p2py)
            self.p2pz = np.array(p2pz)
            self.p2ptot = np.array(p2ptot)
        else:            
            self.p2px = np.array(p2px)[prune_layout[1][:,0]]
            self.p2py = np.array(p2py)[prune_layout[1][:,0]]
            self.p2pz = np.array(p2pz)[prune_layout[1][:,0]]
            self.p2ptot = np.array(p2ptot)[prune_layout[1][:,0]]

        # JobName,Primary,Energy,Zenith,Azimuth,XmaxDistance,SlantXmax,RandomCore[0],RandomCore[1],RandomAzimuth,HadronicModel
        A = showerinfo[0]
        A = A.strip().split()
        self.jobname = A[0]
        self.primary = A[1]
        self.energy = np.float32(A[2])
        self.zenith = 180 - np.float32(A[3])  # CONVERSION TO REAL ZENITH ANGLE
        self.azimuth = np.float32(A[4])
        self.xmax_distance = np.float32(A[5])
        self.slant_xmax = np.float32(A[6])
        self.random_core0 = np.float32(A[7])
        self.random_core1 = np.float32(A[8])
        self.random_azimuth = np.float32(A[9])
        self.hadronic_model = A[10]
        self.step = np.float32(step)
        self.name = name
        self.init_layout()

    def init_layout(self):
        if "hexhex" in self.name:
            self.layout = "hexhex"
        elif "rect" in self.name:
            self.layout = "rect"
        elif "trihex" in self.name:
            self.layout = "trihex"
        else:
            self.layout = "unknown"

    def is_triggered1(self, threshold):
        return self.p2ptot > threshold

def make_ev_list(path):
    ev_list = []
    count = 0

    for subdir in os.listdir(path):
        if os.path.isdir(os.path.join(path, subdir)):
            list_fn = os.listdir(os.path.join(path, subdir)) 
            for fn in list_fn:
                if fn[-6:] == "P2Pdat":
                    f1 = os.path.join(path, subdir, fn)
                    f2 = os.path.join(path, subdir, fn+'.showerinfo')
                    step  = fn.split("_")[-2]

                    ev_list.append(Event(f1, f2, step, fn))

        count += 1 
        if(count % 100 == 0):
            print("Event #{} done".format(count))
    return ev_list

def get_info_from_merged_filename(merged_file):
    bn = os.path.basename(merged_file)
    bn = os.path.splitext(bn)[0]
   
    return bn.split('_')


def make_ev_list_from_merged_file(merged_file, prune_layout=()):
    ev_list = []
    count = 0

    primary, grid_shape, step = get_info_from_merged_filename(merged_file)
    # The following uses ijson as json can't load large files
    with open(merged_file) as f:
        data = ijson.kvitems(f, "", use_float=True, multiple_values=True)
        for k, v in data:
            subdir = k
            if "ef" in v:
                name = subdir + ".Interpolated." + "%s_%s_efield.P2Pdat"%(grid_shape, step) 
                ev_list.append(
                    Event(
                        v["ef_showerinfo"],
                        v["ef"],
                        step, 
                        name,
                        prune_layout
                    )
                )
            
            if "vo" in v: 
                name = subdir + ".Interpolated." + "%s_%s_voltage.P2Pdat"%(grid_shape, step) 
                ev_list.append(
                    Event(
                        v["vo_showerinfo"],
                        v["vo"],
                        step,
                        name,
                        prune_layout
                    )
                )
            
    return ev_list


def make_ev_list_from_merged_file_json(merged_file, prune_layout=()):
    ev_list = []
    count = 0

    primary, grid_shape, step = get_info_from_merged_filename(merged_file)
   
    with open(merged_file) as f:
        data = json.load(f)

    keys = list(data.keys())

    for subdir in keys:
        if "ef" in data[subdir]:
            name = subdir + ".Interpolated." + "%s_%s_efield.P2Pdat"%(grid_shape, step) 
            ev_list.append(
                Event(
                    data[subdir]["ef_showerinfo"],
                    data[subdir]["ef"],
                    step, 
                    name,
                    prune_layout
                )
            )
        if "vo" in data[subdir]:
            name = subdir + ".Interpolated." + "%s_%s_voltage.P2Pdat"%(grid_shape, step) 
            ev_list.append(
                Event(
                    data[subdir]["vo_showerinfo"],
                    data[subdir]["vo"],
                    step,
                    name,
                    prune_layout
                )
            )
            
    return ev_list


def create_ev_select_from_config_merged(
    events_data_dir,
    merged_file_dir,
    config_merged, 
    threshold,
    n_trig_thres,
    prune_layout=()
):
    """
    Creates all the ev_select files for all the configurations defined in the 
    config_merged files, for given threshold and n_trig_thres
    The idea of this function is to have a func that created everything once and for all
    """

    for primary in config_merged["primaries"]:
        for k,v in config_merged["layouts"].items():
            for step in v:
                _ = create_ev_select(
                    events_data_dir,
                    merged_file_dir,
                    k,
                    primary,
                    step,
                    threshold,
                    n_trig_thres,
                    prune_layout
                )


def make_sanity_plots(ev_list, grid_shape, step, primary, sanity_plots_dir, input_n_ring=10):

    pos0, offset0 = grids.create_grid_univ(grid_shape, step, do_prune=False, input_n_ring=input_n_ring)

    offx = []
    offy = []
    angles = []

    for evv in ev_list:
        offx.append(evv.random_core0)
        offy.append(evv.random_core1)
        angles.append(evv.random_azimuth)

    offx = np.array(offx)
    offy = np.array(offy)
    angles = np.array(angles)
    plt.figure() 
    plt.clf()
    plt.hist2d(offx, offy, bins=30)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Distribution of shower core positions')
    plt.colorbar()
    plt.savefig(os.path.join(sanity_plots_dir, "offset_distribution_{}_{}_{}.png".format(grid_shape, step, primary)))

    plt.figure()
    plt.clf()
    plt.hist(angles, bins=30)
    plt.xlabel('Azimuth rotation angle [deg]')
    plt.ylabel('N')
    plt.title('Distribution of random azimuth rotation angles')
    plt.savefig(os.path.join(sanity_plots_dir, "azimuth_distribution_{}_{}_{}.png".format(grid_shape, step, primary)))

    n_ev = len(ev_list)
    dum = np.arange(n_ev)
    np.random.shuffle(dum)
    for i in dum[0:30]:
        k = i
        ev = ev_list[k]

        pos, offset = grids.create_grid_univ(grid_shape, step, do_prune=False, input_n_ring=input_n_ring, angle=ev.random_azimuth)
        off0 = ev.random_core0 
        off1 = ev.random_core1 


        x = off0
        y = off1
        theta = ev.random_azimuth / 180 * np.pi
        xp = x * np.cos(theta) - y * np.sin(theta)
        yp = x * np.sin(theta) + y * np.cos(theta)

        plt.figure()
        plt.xlabel('x [m]')
        plt.ylabel('y [m]') 
        plt.title("%d, used for interp"%k)
        plt.scatter(pos[0]-xp, pos[1]-yp, c=ev.p2ptot)
        plt.axis('equal')
        plt.colorbar()
        plt.savefig(os.path.join(sanity_plots_dir, "rotated_grid_%d.png"%k))
        
        plt.figure()
        plt.xlabel('x [m]')
        plt.ylabel('y [m]') 
        plt.title("%d"%k)
        plt.scatter(pos0[0], pos0[1], c=ev.p2ptot)
        plt.axis('equal')
        plt.colorbar()
        plt.savefig(os.path.join(sanity_plots_dir, "unrotated_grid_%d.png"%k))



def create_ev_select(
    events_data_dir,
    merged_file_dir,
    sanity_plots_dir,
    grid_shape,
    primary,
    step,
    threshold,
    n_trig_thres,
    prune_layout=(), 
    input_n_ring=10
):
    
    # ev_select_name is e.g. ev_select_hexhex_Proton_250_30.0000_5_.npy
    ev_select_name = 'ev_select_%s_%s_%s_%2.4f_%d_%s.npy'
    ev_select_file = os.path.join(events_data_dir, ev_select_name)

    if prune_layout == ():
        ev_select_file = ev_select_file%(
            grid_shape,
            primary,
            step,
            threshold,
            n_trig_thres,
            ""
        )
    else:
        ev_select_file = ev_select_file%(
            grid_shape,
            primary,
            step,
            threshold,
            n_trig_thres,
            prune_layout[0]
        )

    do_make_ev_list = isnot_ev_select(ev_select_file)

    if do_make_ev_list:
        print('creating ev_select_file for {} {} {}'.format(grid_shape, primary, step))
        merged_file = os.path.join(merged_file_dir, '%s_%s_%s.json'%(primary, grid_shape, step))
        ev_list = make_ev_list_from_merged_file(merged_file, prune_layout)
        if prune_layout[0] == "all":
            make_sanity_plots(ev_list, grid_shape, np.float32(step), primary, sanity_plots_dir, input_n_ring=input_n_ring)

        for ev in ev_list:
            if "voltage" in ev.name:
                ev.num_triggered = sum(ev.is_triggered1(threshold))
                ev.is_triggered2 = (ev.num_triggered > n_trig_thres)

        ev_select = [
            (
                ev.num_triggered,
                ev.energy,
                ev.step,
                ev.zenith,
                ev.is_triggered2
            ) for  ev in ev_list
        if "voltage" in ev.name
        ]

        ev_select = np.array(ev_select)
        np.save(ev_select_file, ev_select)
       

def get_ev_select(
    events_data_dir,
    grid_shape,
    primary,
    step,
    threshold,
    n_trig_thres,
    prune_layout=()
):
    
    # ev_select_name is e.g. ev_select_hexhex_Proton_250_30.0000_5_.npy
    ev_select_name = 'ev_select_%s_%s_%s_%2.4f_%d_%s.npy'
    ev_select_file = os.path.join(events_data_dir, ev_select_name)


    if prune_layout == ():
        ev_select_file = ev_select_file%(
            grid_shape,
            primary,
            step,
            threshold,
            n_trig_thres,
            ""
        )
    else:
        ev_select_file = ev_select_file%(
            grid_shape,
            primary,
            step,
            threshold,
            n_trig_thres,
            prune_layout[0]
        )

   
    ev_select = np.load(ev_select_file)
    return ev_select


def isnot_ev_select(ev_select_file):
    return not(os.path.isfile(ev_select_file))


def make_ev_select_old(ev_list, layout, primary, ev_select_file):
    ev_select = [
        (
            ev.num_triggered,
            ev.energy,
            ev.step,
            ev.zenith,
            ev.is_triggered2
        ) for  ev in ev_list
        if "voltage" in ev.name
        and ev.primary == primary
        and ev.layout == layout
    ]
    ev_select = np.array(ev_select)
    np.save(ev_select_file, ev_select)
    return ev_select


def compute_meanNtrig(enerbins_limits, zenbins_limits, ev_select):
    """
    compute the mean of variance of number of triggered antennas
    in various energy and zenith bins

    """
   
    n_energy_bins = len(enerbins_limits) - 1
    n_zenith_bins = len(zenbins_limits) - 1

    mean_n_trig = np.zeros((n_energy_bins, n_zenith_bins))
    var_n_trig = np.zeros((n_energy_bins, n_zenith_bins))
   
    meanNtrig_ener = []
    varNtrig_ener = []


    for i_ener in range(n_energy_bins):
        for i_zen in range(n_zenith_bins):
            ind = np.where(
                (ev_select[:,1] > enerbins_limits[i_ener]) * 
                (ev_select[:,1] <= enerbins_limits[i_ener+1]) *
                (ev_select[:,3] > zenbins_limits[i_zen]) *
                (ev_select[:,3] <= zenbins_limits[i_zen+1])
            )

            mean_n_trig[i_ener, i_zen] = np.mean(ev_select[ind[0],0])
            var_n_trig[i_ener, i_zen] = np.var(ev_select[ind[0],0])
    
    return mean_n_trig, var_n_trig


def compute_trig_rate(enerbins_limits, zenbins_limits, ev_select):
    """
    compute the number of triggering events in various energy and zenith bins
    """
    n_energy_bins = len(enerbins_limits) - 1
    n_zenith_bins = len(zenbins_limits) - 1

    rate_trig = np.zeros((n_energy_bins, n_zenith_bins))


    for i_ener in range(n_energy_bins):
        for i_zen in range(n_zenith_bins):
            ind = np.where(
                (ev_select[:,1] > enerbins_limits[i_ener]) * 
                (ev_select[:,1] <= enerbins_limits[i_ener+1]) *
                (ev_select[:,3] > zenbins_limits[i_zen]) *
                (ev_select[:,3] <= zenbins_limits[i_zen+1])
            )
            if len(ind[0]) == 0:
                rate_trig[i_ener, i_zen] = 0
            else:
                rate_trig[i_ener, i_zen] = sum(ev_select[ind[0],4])/np.size(ev_select[ind[0],4])
    return rate_trig

    # for iener, ener in enumerate(enerbins):
    #         Ntrig2_zen = []
        
    #         #for izen in range(0, len(zenbins)-1):
    #         for izen, zen in enumerate(zenbins):
    #             ind = np.where((ev_select[:,1] == ener) * (ev_select[:,2] == step) 
    #                 *(np.abs(ev_select[:,3]-(180-zen)) < 0.5))
    
    #             if len(ind[0])==0:
    #                 Ntrig2_zen.append(0)
    #             else:
    #                 Ntrig2_zen.append(sum(ev_select[ind[0],4])/np.size(ev_select[ind[0],4]))

    #         Ntrig2_step.append(Ntrig2_zen)

    #     Ntrig2_ener.append(Ntrig2_step)

    # Ntrig2_ener = np.array(Ntrig2_ener)
    # return Ntrig2_ener


def plot_meanNtrig(
    meanNtrig_ener,
    varNtrig_ener,
    enerbins,
    enerbins_center,
    zenbins,
    zenbins_center,
    layout="trihex"
):
    """
    plot Ntriggered antennas vs energies for fixed steps
    enerbins and zenbins are the arrays of the bin limits
    """
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




def plot_trigrate(
    trig_rate,
    enerbins,
    enerbins_center,
    zenbins,
    zenbins_center,
    layout="trihex"
):
    """
    plot Ntriggered antennas vs energies for fixed steps
    enerbins and zenbins are the arrays of the bin limits
    """
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



# def plot_Ntrig_fixedzenith_vsenergy(
#     meanNtrig_ener,
#     varNtrig_ener,
#     stepbins,
#     enerbins,
#     zenbins,
#     layout="rect"
# ):
#     """
#     plot Ntriggered antennas vs energies for fixed zenith angles
#     """
#     for izen in range(0, len(zenbins)-1):
#         plt.figure(izen+4) 
#         plt.clf()
#         for istep, step in enumerate(stepbins):
#             plt.errorbar(
#                 enerbins,
#                 meanNtrig_ener[istep,:,izen],
#                 yerr=np.sqrt(varNtrig_ener[istep,:,izen]), 
#                 fmt=SYM_LIST[istep],
#                 capsize=2,
#                 alpha=0.7,
#                 label='step = %d m'%(np.int32(step)))
#             #plt.errorbar(enerbins, Ntrig2_ener[istep,:,izen],  
#             #   fmt=sym_list[istep], capsize=2, alpha=0.7)
#         plt.yscale('log')
#         plt.ylabel('N triggered antennas')
#         plt.xlabel('energy [EeV]')
#         plt.title('%s, %4.0f > zenith >%4.0f deg'%(layout, 180-zenbins[izen], 180-zenbins[izen+1]))
#         plt.legend(loc=4)
#         #plt.show()


# def plot_Ntrig_fixedernergy_vszenith(
#     meanNtrig_ener,
#     varNtrig_ener,
#     stepbins,
#     enerbins,
#     zenbins,
#     layout="rect",
#     plot_path='./', 
#     primary="Proton"
# ):
#     """
#     # plot Ntriggered antennas vs zenith angles for fixed energies 
#     """
#     for iener, ener in enumerate(enerbins):
#         plt.figure(iener) 
#         plt.clf()
#         for istep, step in enumerate(stepbins):
#             plt.errorbar(
#                 zenbins,
#                 meanNtrig_ener[istep,iener,:],
#                 yerr=np.sqrt(varNtrig_ener[istep,iener,:]), 
#                 fmt=SYM_LIST[istep],
#                 capsize=2,
#                 alpha=0.7,
#                 label='step = %d m'%(np.int32(step))
#             )
#             #plt.errorbar(enerbins, Ntrig2_ener[istep,:,izen], 
#             #    fmt=sym_list[izen], capsize=2)
#         plt.yscale('log')
#         plt.ylabel('N triggered antennas')
#         plt.xlabel('zenith [deg]')
#         plt.title('%s, %s, E = %4.3f EeV'%(primary,layout, ener))
#         plt.legend(loc=2)
#         plt.ylim(1,225)
#         plt.xlim(45,90)
#         #plt.show()
#         plt.savefig(os.path.join(plot_path, 'Ntrig_vs_zen_E%4.3f_%s_%s.png'%(ener, layout, primary)))


# def plot_rate_fixedsteps_vsenergy(
#     Ntrig2_ener,
#     stepbins,
#     enerbins,
#     zenbins,
#     layout="rect",
#     plot_path='./'
# ):
#     for istep, step in enumerate(stepbins):
#         plt.figure(istep) 
#         plt.clf()
#         for izen in range(0, len(zenbins)-1):
#             plt.errorbar(
#                 enerbins,
#                 Ntrig2_ener[istep,:,izen], 
#                 fmt=SYM_LIST[izen],
#                 ls='-',
#                 capsize=2,
#                 alpha=0.7,
#                 label='%4.0f > zen >%4.0f deg'%(180-zenbins[izen], 180-zenbins[izen+1])
#             )
#         plt.yscale('log')
#         plt.ylabel('Triggered event rate')
#         plt.xlabel('energy [EeV]')
#         plt.title('%s, step = %d m'%(layout, np.int32(step)))
#         plt.legend(loc=4)
#         #plt.show()
#         plt.savefig(os.path.join(plot_path,'trigevrate_vs_energy_step%d_%s_30muV.png'%(np.int32(step), layout)))


# def plot_rate_fixedzenith_vsenergy(
#     Ntrig2_ener,
#     stepbins,
#     enerbins,
#     zenbins,
#     layout="rect",
#     plot_path="./"   
# ):
#     for izen in range(0, len(zenbins)-1):
#         plt.figure(izen+4) 
#         plt.clf()
#         for istep, step in enumerate(stepbins):
#             plt.errorbar(
#                 enerbins,
#                 Ntrig2_ener[istep,:,izen],  
#                 fmt=SYM_LIST[istep],
#                 ls='-',
#                 capsize=2,
#                 alpha=0.7,
#                 label='step = %d m'%(np.int32(step))
#             )
#         plt.yscale('log')
#         plt.ylabel('Triggered event rate')
#         plt.xlabel('energy [EeV]')
#         plt.title('%s, %4.0f > zenith >%4.0f deg'%(layout, 180-zenbins[izen], 180-zenbins[izen+1]))
#         plt.legend(loc=4)
#         #plt.show()
#         plt.savefig(os.path.join(plot_path,'trigevrate_vs_energy_z%4.1f_%s_30muV.png'%(180-zenbins[izen+1], layout)))


# def plot_rate_fixedenergy_vszenith(
#     Ntrig2_ener,
#     stepbins,
#     enerbins,
#     zenbins,
#     layout="rect",
#     plot_path="./"  
# ):
# # plot Ntriggered events vs zenith angles for fixed energies 
#     for iener, ener in enumerate(enerbins):
#         plt.figure(iener) 
#         plt.clf()
#         for istep, step in enumerate(stepbins):
#             plt.errorbar(
#                 zenbins,
#                 Ntrig2_ener[istep,iener,:],  
#                 fmt=SYM_LIST[istep],
#                 ls='-',
#                 capsize=2,
#                 alpha=0.7,
#                 label='step = %d m'%(np.int32(step))
#             )
#         plt.yscale('log')
#         plt.ylabel('Triggered event rate')
#         plt.xlabel('zenith [deg]')
#         plt.title('Proton, %s, E = %4.3f EeV'%(layout, ener))
#         plt.legend(loc=4)
#         plt.ylim(1.e-2,1.1)
#         plt.xlim(45,90)
#         #plt.show()
#         plt.savefig(os.path.join(plot_path, 'trigevrate_vs_zen_E%4.3f_%s_Proton.png'%(ener, layout)))
