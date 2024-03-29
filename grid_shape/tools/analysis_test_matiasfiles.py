import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 



path = "/Users/benoitl/Documents/GRAND/GridsTest/"

#plot_path = '/Users/kotera/BROQUE/Plots_GRAND/'
plot_path = '/Users/benoitl/Documents/GRAND/plots_GridsTest'
#events_data_dir = "/Users/kotera/BROQUE/Data_GRAND/Matias/event_data"
events_data_dir = "/Users/benoitl/Documents/GRAND/event_data_GridsTest"

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



grid_shape = "hexhex"
radius = 125



pos, offset = grids.create_grid_univ(grid_shape, radius, do_prune=False)


trihex_steps = config_merged["layouts"]["trihex"]
hexhex_steps = config_merged["layouts"]["hexhex"]


[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        grid_shape,
        "Gamma",
        step, 
        threshold,
        n_trig_thres
    )
    for step in hexhex_steps
]



ev_select_hexhex = [
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

ev_select_hexhex = np.concatenate([*ev_select_hexhex])  


ev_select_trihex = [
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
ev_select_trihex = np.concatenate([*ev_select_trihex])  



ev_select = ev_select_trihex
layout = "trihex"

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
enerbins = np.unique(ev_select[:,1])
#zenbins = 180-np.unique(A_rect[:,3])
zenbins = np.array([94.77,95.74,97.18,98.21,99.59,101.54, 104.48, 106.6, 109.47, 113.58,120,132])
zenbins = 180. - zenbins
#zenbins = [94,100,105,110,120,131]
stepbins = np.unique(ev_select[:,2])

meanNtrig_ener1, varNtrig_ener1 = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_trihex
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

Ntrig2_ener_hexhex1 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_trihex
)


rate_hexhex1 = Ntrig2_ener_hexhex1.copy() * 0

rate_hexhex1_area = Ntrig2_ener_hexhex1.copy() * 0


for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            
            rate_hexhex1[istep, iener, izen] = (
                Ntrig2_ener_hexhex1[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_trihex[istep] * np.cos(zen*np.pi/180) /area_trihex[istep]*1e6 # !! finally we divide by the area here
            )
    
            rate_hexhex1_area[istep, iener, izen] = (
                Ntrig2_ener_hexhex1[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_trihex[istep] * np.cos(zen*np.pi/180) # 
            )
            



for izen in range(0, len(zenbins)-1):
    plt.figure(izen) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        
        plt.errorbar(
            enerbins,
            rate_hexhex1[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'h',
            ms=10,           
            ls='-',
            capsize=2,
            alpha=0.7,
            #label='step = %d m'%(np.int32(step))
        )
       

    plt.yscale('log')
    plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} km^{-2}]$')
    plt.xlabel('energy [EeV]')
    plt.title('%s, %4.2f > zenith > %4.2f deg'%("trihex", thetar[izen], thetal[izen+1]))
    plt.legend(loc=1)
    #plt.ylim(1.e-1,1.e2)
    #plt.ylim(1.e-5,1.e1)        
    #plt.show()
    #plt.savefig('compmethod_evrate_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1]))



for izen in range(0, len(zenbins)-1):
    plt.figure(izen) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        
        plt.errorbar(
            enerbins,
            rate_hexhex1_area[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'h',
            ms=10,           
            ls='-',
            capsize=2,
            alpha=0.7,
            #label='step = %d m'%(np.int32(step))
        )
       

    plt.yscale('log')
    plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1}]$')
    plt.xlabel('energy [EeV]')
    plt.title('%s, %4.2f > zenith > %4.2f deg'%("trihex", thetar[izen], thetal[izen+1]))
    plt.legend(loc=1)
    #plt.ylim(1.e-1,1.e2)
  #  plt.ylim(1.e-4,1.e2)        
    #plt.show()
   # plt.savefig('compmethod_area_evrate_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1]))






plt.figure(1)
plt.clf()
plt.scatter(pos[0,:], pos[1,:], c = mask[:,0])
plt.axis("equal")
plt.colorbar()






plt.show()






file_trihex = "/Users/benoitl/Documents/GRAND/StshpLibrary-HDF5-Grids_merge/Proton_trihex_125.json"
file_hexhex = "/Users/benoitl/Documents/GRAND/StshpLibrary-HDF5-Grids_merge/Proton_hexhex_125.json"

with open(file_trihex) as f:
    trihex_data = json.load(f)


with open(file_hexhex) as f:
    hexhex_data = json.load(f)


dirr = "Stshp_XmaxLibrary_0.03162_82.82_90_Proton_10"
dirr = "Stshp_XmaxLibrary_0.03981_73.39_90_Proton_05"
hx,hy,hz,ht = hexhex_data[dirr]['ef']

tx,ty,tz,tt = trihex_data[dirr]['ef']



pos_hex, _ = grids.create_grid_univ("hexhex", radius)
pos_tri, _ = grids.create_grid_univ("trihex", radius)



plt.figure(3)
plt.clf()
plt.scatter(pos_hex[0], pos_hex[1], c = ht)
plt.axis("equal")
plt.colorbar()

plt.figure(4)
plt.clf()
plt.scatter(pos_tri[0], pos_tri[1], c = tt)
plt.axis("equal")
plt.colorbar()






