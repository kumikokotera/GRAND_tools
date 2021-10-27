import numpy as np
import os
import json
import matplotlib.pyplot as plt
from interpolation_functions import hdf5fileinout as hd
import ijson
from grid_shape_lib import grids as grids

'''
Definition of class to read the simulation output files
and of event selection routines
to perform simulation output analysis in the analysis.py routine
'''

SYM_LIST = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
MYC = ['0','0.20','0.4','0.6','0.8']


class Event_old:
    def __init__(self, f1, f2, step, name):
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
            self.p2px = np.array(p2px)[prune_layout[1][:, 0]]
            self.p2py = np.array(p2py)[prune_layout[1][:, 0]]
            self.p2pz = np.array(p2pz)[prune_layout[1][:, 0]]
            self.p2ptot = np.array(p2ptot)[prune_layout[1][:, 0]]
        
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
    """ Deprecated? """

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


def make_ev_list_from_hdf5_files(events_data_dir):
    files_list = []
    ev_list = []
    count = 0

    # list all files 
    for path, sub, files in os.walk(events_data_dir):
        [files_list.append(os.path.join(path, f)) for f in files if os.path.splitext(f)[-1] == '.hdf5' ]
    import time
    t0 = time.time()
    for f in files_list:
        count += 1
        print(count)
        run_info = hd.GetRunInfo(f)
        event_name = run_info["EventName"][0]
        #print(event_name)
        event_info = hd.GetEventInfo(f, event_name)
        #print(event_info)
        rand_azimuth = os.path.splitext(os.path.basename(f))[0].split('_')[-1]
        # JobName,Primary,Energy,Zenith,Azimuth,XmaxDistance,SlantXmax,RandomCore[0],RandomCore[1],RandomAzimuth,HadronicModel
        info_string = "{},{},{},{},{},{},{},{},{},{},{}".format(
            event_name,
            event_info["Primary"],
            event_info["Energy"],
            event_info["Zenith"],
            event_info["Azimuth"],
            event_info["XmaxDistance"],
            event_info["SlantXmax"],
            event_info["CorePosition"][0][0],
            event_info["CorePosition"][0][1],
            rand_azimuth,
            run_info["HadronicModel"][0]
        )
        antenna_info = hd.GetAntennaP2PInfo(f, event_name[0])
        antenna_id = antenna_info["ID"]
        antenna_p2px = antenna_info[""]

        #ßshower_data = 

    print(time.time() - t0)
    return files_list


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
    """
    Deprecated?
    """
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
    energies = []
    zeniths = []
    
    for evv in ev_list:
        offx.append(evv.random_core0)
        offy.append(evv.random_core1)
        angles.append(evv.random_azimuth)
        energies.append(evv.energy)
        zeniths.append(evv.zenith)
        

    offx = np.array(offx)
    offy = np.array(offy)
    angles = np.array(angles)
    energies = np.array(energies)
    zeniths  = np.array(energies)
    

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
        plt.title(
            "%2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f"%(
                k,
                ev.random_core0,
                ev.random_core1,
                ev.random_azimuth,
                ev.zenith,
                ev.energy
            ),
            fontdict={"fontsize":10}
        )
        plt.scatter(pos[0]-xp, pos[1]-yp, c=ev.p2ptot)
        plt.axis('equal')
        plt.colorbar()
        plt.savefig(os.path.join(sanity_plots_dir, "rotated_grid_%d.png"%k))
        
        plt.figure()
        plt.xlabel('x [m]')
        plt.ylabel('y [m]') 
        #plt.title("%d"%k)
        plt.title(
            "%2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f"%(
                k,
                ev.random_core0,
                ev.random_core1,
                ev.random_azimuth,
                ev.zenith,
                ev.energy
            ),
            fontdict={"fontsize":10}
        )
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
    input_n_ring=10, 
    hdf5=False
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
        if not(hdf5):
            merged_file = os.path.join(merged_file_dir, '%s_%s_%s.json'%(primary, grid_shape, step))
            ev_list = make_ev_list_from_merged_file(merged_file, prune_layout)
        else:
            ev_list = make_ev_list_from_hdf5_files()
        for ev in ev_list:
            if "voltage" in ev.name:
                ev.num_triggered = sum(ev.is_triggered1(threshold))
                ev.is_triggered2 = (ev.num_triggered > n_trig_thres)
        
        if prune_layout[0] == "all":
            make_sanity_plots(ev_list, grid_shape, np.float32(step), primary, sanity_plots_dir, input_n_ring=input_n_ring)

        ev_select = [
            (
                ev.num_triggered,
                ev.energy,
                ev.step,
                ev.zenith,
                ev.is_triggered2,
                ev.azimuth,
                ev.random_core0, ##add
                ev.random_core1, ##add
                ev.random_azimuth ##add
            ) for ev in ev_list
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
