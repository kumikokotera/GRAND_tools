import numpy as np
import matplotlib.pyplot as plt
import os
import json

from grid_shape_lib.modules import grids as grids
from grid_shape_lib.modules import hexy as hx

from grid_shape_lib.utils import utils_analysis as ua
from grid_shape_lib.utils import diff_spec


Z_SITE = 2900 # height of GP300 site in km


## contruct all hexes up to n_ring = 4 and center at n=0, none at n=1, 
# one out of two at n=2, none at n=3, and one out of two at n=4

radius = 1000


#hex4, _ = grids.create_grid_univ('hexhex', radius, input_n_ring=4)
def get_triisle_mask(n_ring, radius, do_plots=False):
 
    hexx, _ = grids.get_hexarray(n_ring, radius)
    
    n = 0
    xcube0 = hx.get_disk((0,0,0), n)
    xpix0 = hx.cube_to_pixel(xcube0, radius)

    n = 2 
    xcube2 = hx.get_ring((0,0,0), n)[::2,:]
    xpix2 = hx.cube_to_pixel(xcube2, radius)

    n = 4
    xcube4 = hx.get_ring((0,0,0), n)[::2,:]
    xpix4 = hx.cube_to_pixel(xcube4, radius)

    pos = np.vstack([hexx, xpix0, xpix2, xpix4])

    if do_plots:
        plt.figure(1)
        plt.clf()
        plt.plot(hexx[:,0], hexx[:,1], 'k.')
        plt.plot(xpix0[:,0], xpix0[:,1], 'r.')
        plt.plot(xpix2[:,0], xpix2[:,1], 'g.')
        plt.plot(xpix4[:,0], xpix4[:,1], 'c.')
        plt.plot(pos[:,0]+100, pos[:,1]+100, 'm.')
        plt.axis('equal')

    grid_x = pos[:,0]
    grid_y = pos[:,1]

        # remove redundant points
    x_pos_new, y_pos_new = grids.remove_redundant_point(grid_x, grid_y)
    n_pos_new = len(x_pos_new)

    tri, _ = grids.create_grid_univ('trihex', radius, input_n_ring=n_ring)
    n_tri = len(tri[0])
    mask = np.zeros((n_tri, 1))

    t = [(np.int32(tri[0,i]*10000), np.int32(tri[1,i]*10000)) for i in range(n_tri)]
    tnew = [(np.int32(x_pos_new[i]*10000), np.int32(y_pos_new[i]*10000)) for i in range(n_pos_new)]

    for i in range(n_tri):
        if t[i] in tnew:
            mask[i] = 1 

    mask = mask.astype(bool)
    if do_plots:
        plt.figure(2, figsize=(8, 5))
        plt.clf()
        plt.scatter(tri[0], tri[1], c = 1-mask[:,0])
        a = plt.axis("equal")
        plt.title('Triisle4')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]') 
        plt.savefig("triisle4.png")
    return mask


## faire trihex4
tri4, _ = grids.create_grid_univ("trihex", 1000, input_n_ring=4)
hex4, _, m4 = grids.create_grid_univ("trihex", 1000, input_n_ring=4, do_prune=True)

plt.figure(3, figsize=(8, 5))
plt.clf()
plt.scatter(tri4[0], tri4[1], c = 0*m4[:,0]+1)
a = plt.axis("equal")
plt.title('trihex4')
plt.xlabel('x [m]')
plt.ylabel('y [m]') 
plt.savefig("trihex4.png")

plt.figure(4, figsize=(8, 5))
plt.clf()
plt.scatter(hex4[0], hex4[1], c = 1-m4[:,0])
a = plt.axis("equal")
plt.title('hexhex4')
plt.xlabel('x [m]')
plt.ylabel('y [m]') 
plt.savefig("hexhex4.png")









# hexarray, radius_grid,  mask = get_hexarray(n_ring, radius, do_mask=True)

# grid_x = hexarray[:,0]
# grid_y = hexarray[:,1]

#         # remove redundant points
# x_pos_new, y_pos_new, mask = remove_redundant_point(grid_x, grid_y, mask=mask)
    





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



# make the masks
_, _, mask1 = grids.create_grid_univ("trihex", radius, do_prune=True, input_n_ring = 4)
mask2 = get_triisle_mask(4, 1000)


trihex_steps = config_merged["layouts"]["trihex"]
# create the ev_selects for the masks
[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        "trihex",
        "Proton",
        step, 
        threshold,
        n_trig_thres,
        prune_layout=("t_to_h4", mask1)
    )
    for step in trihex_steps
]

[
    ua.create_ev_select(
        events_data_dir,
        merged_file_dir,
        "trihex",
        "Proton",
        step, 
        threshold,
        n_trig_thres,
        prune_layout=("t_to_triisle4", mask2)
    )
    for step in trihex_steps
]

## 1 faire les plots trihex4 purs
## 2 faire les plots trihex4 pruned to h4
## 3 faire les plots thihex4 pruned to 36326263-4


ev_select_t_to_h4= [
    ua.get_ev_select(
        events_data_dir,
        "trihex",
        "Proton",
        step,
        threshold,
        n_trig_thres,
        prune_layout=("t_to_h4", mask1)
    )
    for step in trihex_steps
]
ev_select_t_to_h4 = np.concatenate([*ev_select_t_to_h4])  

ev_select_t_to_triisle4 = [
    ua.get_ev_select(
        events_data_dir,
        "trihex",
        "Proton",
        step,
        threshold,
        n_trig_thres,
        prune_layout=("t_to_triisle4", mask2)
    )
    for step in trihex_steps
]
ev_select_t_to_triisle4 = np.concatenate([*ev_select_t_to_triisle4])  


ev_select_t4 = [
    ua.get_ev_select(
        events_data_dir,
        "trihex",
        "Proton",
        step,
        threshold,
        n_trig_thres,
        prune_layout=()
    )
    for step in trihex_steps
]
ev_select_t4 = np.concatenate([*ev_select_t4])  




ev_select = ev_select_t_to_h4
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
    ev_select_t_to_h4
)



meanNtrig_ener2, varNtrig_ener2 = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_t_to_triisle4
)



meanNtrig_ener3, varNtrig_ener3 = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_t4
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

Ntrig2_ener_t_to_h4 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_t_to_h4
)


Ntrig2_ener_t_to_triisle4 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_t_to_triisle4
)
Ntrig2_ener_t4 = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_t4
)

rate_t_to_h4 = Ntrig2_ener_t_to_h4.copy() * 0
rate_t_to_triisle4 = Ntrig2_ener_t_to_triisle4.copy() * 0
rate_t4 = Ntrig2_ener_t4.copy() * 0

rate_t_to_h4_area = Ntrig2_ener_t_to_h4.copy() * 0
rate_t_to_triisle4_area = Ntrig2_ener_t_to_triisle4.copy() * 0
rate_t4_area = Ntrig2_ener_t4.copy() * 0



for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            
            rate_t_to_h4[istep, iener, izen] = (
                Ntrig2_ener_t_to_h4[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) /area_trihex[istep]*1e6 # !! finally we divide by the area here
            )
            rate_t_to_triisle4[istep, iener, izen] = (
                Ntrig2_ener_t_to_triisle4[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) /area_trihex[istep]*1e6 # !! finally we divide by the area here
            )
            rate_t4[istep, iener, izen] = (
                Ntrig2_ener_t4[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) /area_trihex[istep]*1e6 # !! finally we divide by the area here
            )
            rate_t_to_h4_area[istep, iener, izen] = (
                Ntrig2_ener_t_to_h4[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) 
            )
            rate_t_to_triisle4_area[istep, iener, izen] = (
                Ntrig2_ener_t_to_triisle4[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) 
            )
            rate_t4_area[istep, iener, izen] = (
                Ntrig2_ener_t4[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) 
            )



for izen in range(0, len(zenbins)-1):
    plt.figure(izen) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        
        plt.errorbar(
            enerbins,
            rate_t_to_h4[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'h',
            ms=10,           
            ls='-',
            capsize=2,
            alpha=0.7,
            #label='step = %d m'%(np.int32(step))
        )
        plt.errorbar(
            enerbins,
            rate_t_to_triisle4[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'.',
            ms=10,           
            ls=':',
            capsize=2,
            alpha=0.7,
            #label='step = %d m'%(np.int32(step))
        )
        plt.errorbar(
            enerbins,
            rate_t4[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'^',
            ms=8,           
            ls='--',
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
    plt.savefig('evrate_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1]))




for izen in range(0, len(zenbins)-1):
    plt.figure(izen) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        
        plt.errorbar(
            enerbins,
            rate_t_to_h4_area[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'h',
            ms=10,           
            ls='-',
            capsize=2,
            alpha=0.7,
            #label='step = %d m'%(np.int32(step))
        )
        plt.errorbar(
            enerbins,
            rate_t_to_triisle4_area[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'.',
            ms=10,           
            ls=':',
            capsize=2,
            alpha=0.7,
            #label='step = %d m'%(np.int32(step))
        )
        plt.errorbar(
            enerbins,
            rate_t4_area[istep,:,izen] * 24*3600,
            fmt="C%d"%istep+'^',
            ms=8,           
            ls='--',
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
    plt.ylim(1.e-5, 1.e2)        
    #plt.show()
    plt.savefig('evrate_area_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1]))

