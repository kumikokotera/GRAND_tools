import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 




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



pos, offset, mask = grids.create_grid_univ(grid_shape, radius, do_prune=True)


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





##### try to remove the outer ring in hexhex

hex4, _ = grids.create_grid_univ('hexhex', 1000, input_n_ring=4)
hex5, _ = grids.create_grid_univ('hexhex', 1000, input_n_ring=5)

mask5 = np.zeros((216, 1))


h4 = [(np.int32(hex4[0,i]*10000), np.int32(hex4[1,i]*10000)) for i in range(150)]
h5 = [(np.int32(hex5[0,i]*10000), np.int32(hex5[1,i]*10000)) for i in range(216)]

for i in range(216):
    if h5[i] in h4:
        mask5[i] = 1 

mask5 = mask5.astype(bool)

plt.figure(245)
plt.clf()
plt.scatter(hex5[0], hex5[1], c = mask5[:,0])




[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        "hexhex",
        "Proton",
        step, 
        threshold,
        n_trig_thres,
        prune_layout=("prune5to4", mask5)
    )
    for step in hexhex_steps
]



ev_select_hexhex5to4 = [
    ua.get_ev_select(
        events_data_dir,
        "hexhex",
        "Proton",
        step,
        threshold,
        n_trig_thres,
        prune_layout=("prune5to4", mask5)
    )
    for step in hexhex_steps
]
ev_select_hexhex5to4 = np.concatenate([*ev_select_hexhex5to4])  



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



meanNtrig_ener2, varNtrig_ener2 = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_hexhex5to4
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


Ntrig2_ener_hexhex2 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_hexhex5to4
)

rate_hexhex1 = Ntrig2_ener_hexhex1.copy() * 0
rate_hexhex2 = Ntrig2_ener_hexhex2.copy() * 0

rate_hexhex1_area = Ntrig2_ener_hexhex1.copy() * 0
rate_hexhex2_area = Ntrig2_ener_hexhex2.copy() * 0


for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            
            rate_hexhex1[istep, iener, izen] = (
                Ntrig2_ener_hexhex1[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_trihex[istep] * np.cos(zen*np.pi/180) /area_trihex[istep]*1e6 # !! finally we divide by the area here
            )
            rate_hexhex2[istep, iener, izen] = (
                Ntrig2_ener_hexhex2[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) /area_hexhex[istep]*1e6# !! finally we divide by the area here
            )
            rate_hexhex1_area[istep, iener, izen] = (
                Ntrig2_ener_hexhex1[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_trihex[istep] * np.cos(zen*np.pi/180) # 
            )
            rate_hexhex2_area[istep, iener, izen] = (
                Ntrig2_ener_hexhex2[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) # 
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
        plt.errorbar(
            enerbins,
            rate_hexhex2[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'X',
            ms=10,           
            ls=':',
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
    plt.ylim(1.e-5,1.e1)        
    #plt.show()
    plt.savefig('compmethod_evrate_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1]))



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
        plt.errorbar(
            enerbins,
            rate_hexhex2_area[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'X',
            ms=10,           
            ls=':',
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
    plt.ylim(1.e-4,1.e2)        
    #plt.show()
    plt.savefig('compmethod_area_evrate_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1]))






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






