import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 



path = "/Users/benoitl/Documents/GRAND/test_pruning/"

#plot_path = '/Users/kotera/BROQUE/Plots_GRAND/'
plot_path = '/Users/benoitl/Documents/GRAND/test_pruning/plots'
#events_data_dir = "/Users/kotera/BROQUE/Data_GRAND/Matias/event_data"
events_data_dir = "/Users/benoitl/Documents/GRAND/test_pruning"

os.makedirs(plot_path, exist_ok=True)
os.makedirs(events_data_dir, exist_ok=True)

threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 

merged_file_dir = os.path.join(path, 'merged_dir')

config_json_file = os.path.join(merged_file_dir, 'merge_config.json')

with open(config_json_file) as f:
    config_merged = json.load(f)










pos, offset = grids.create_grid_univ("hexhex", 125, do_prune=False, input_n_ring=4)
# Make pruning mask
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 125, do_prune=True, input_n_ring=4)
pos3, offset3 = grids.create_grid_univ("trihex", 125, do_prune=False, input_n_ring=4)


plt.figure()
plt.scatter(pos2[0], pos2[1], c = mask2[:,0], s=1)


trihex_steps = config_merged["layouts"]["trihex"]
hexhex_steps = config_merged["layouts"]["hexhex"]
grid_shape = "trihex"

## Create ev_select with pruning 
[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        "trihex",
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
        "hexhex",
        "Gamma",
        step, 
        threshold,
        n_trig_thres,
    )
    for step in hexhex_steps
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



ev_select_hexhex_nopruning = [
    ua.get_ev_select(
        events_data_dir,
        "hexhex",
        "Gamma",
        step,
        threshold,
        n_trig_thres
    )
    for step in hexhex_steps
]
ev_select_hexhex_nopruning= np.concatenate([*ev_select_hexhex_nopruning])  



plt.figure()
plt.hist(ev_select_trihex_pruning[:,0], bins = 50, range=[0, 110], label = 'No pruning')
plt.hist(ev_select_hexhex_nopruning[:,0], bins = 50, range=[0, 110], alpha=0.5, label = 'pruning')
plt.ylabel("# of event")
plt.xlabel('Num_triggered')
plt.legend(loc=0)





ev_select = ev_select_hexhex_nopruning
layout = "trihex"

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
enerbins = np.unique(ev_select[:,1])
enerbins = np.array([enerbins[1]])
zenbins = 180-np.unique(ev_select[:,3])
#zenbins = np.array([94.77,95.74,97.18,98.21,99.59,101.54, 104.48, 106.6, 109.47, 113.58,120,132])
#zenbins = np.array([94.77, 132])
#zenbins = 180. - zenbins
#zenbins = [94,100,105,110,120,131]
stepbins = np.unique(ev_select[:,2])

meanNtrig_ener1, varNtrig_ener1 = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_trihex_pruning
)


meanNtrig_ener2, varNtrig_ener2 = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_hexhex_nopruning
)






# trigger rate calculation over full array
# convolving with measured CR flux

log_E_eV = np.log10(enerbins*1.e18) -0.05
log_E_eV = np.append(log_E_eV, log_E_eV[-1]+0.05)
enerbins2 = 10**(log_E_eV) / 1e18 # in EeV
delta_E = enerbins2[1:] - enerbins2[:-1]


cen = 1.0/np.cos(zenbins*np.pi/180)
cenl = cen + 0.25
cenr = cen - 0.25

thetal = np.arccos(1.0/cenl) * 180/np.pi
thetar = np.arccos(1.0/cenr) * 180/np.pi


# delta_omega = - (zenbins[1:] - zenbins[:-1])
# delta_omega = np.insert(delta_omega, 0, delta_omega[0])

delta_theta = thetal - thetar
delta_omega = 2*np.pi * delta_theta *np.pi/180 * np.sin(np.pi/2 - zenbins*np.pi/180)

area_hexhex = stepbins**2 *3*np.sqrt(3)/2 * 91
area_trihex = stepbins**2 *np.sqrt(3)/4 * 61 * 6
# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_1 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_trihex_pruning
)


Ntrig2_ener_2 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_hexhex_nopruning
)


ua.plot_Ntrig_fixedernergy_vszenith(
    meanNtrig_ener2,
    varNtrig_ener2,
    stepbins,
    enerbins,
    zenbins,
    layout="hexhex_nopruning",
    plot_path=plot_path,
    primary="Gamma"
)


ua.plot_Ntrig_fixedernergy_vszenith(
    meanNtrig_ener1,
    varNtrig_ener1,
    stepbins,
    enerbins,
    zenbins,
    layout="trihex_pruning",
    plot_path=plot_path,
    primary="Gamma"
)