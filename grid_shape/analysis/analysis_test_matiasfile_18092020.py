import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
from grid_shape import layout as layout

### Commented analysis script

path = "/Users/benoitl/Documents/GRAND/Data_grids/20200918/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/Trihex"


threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event
 
 
### "Creating" the layout
pos, offset = grids.create_grid_univ("trihex", 125, do_prune=False, input_n_ring=10)

## Creating the pruned layout
pos2, offset2, mask2 = grids.create_grid_univ("trihex", 125, do_prune=True, input_n_ring=10)

mask  = mask2 + True

lay1 = layout.Layout(path, pos, mask, "all", threshold, n_trig_thres)
lay2 = layout.Layout(path, pos2, mask2, "simple", threshold, n_trig_thres)





##  plot showing the histogram of number of triggered antenna accross the events.
plt.figure()
plt.hist(lay1.ev_select[:,0], bins = 500, range=[0, 1200], label = 'No pruning')
plt.hist(lay2.ev_select[:,0], bins = 500, range=[0, 1200], alpha=0.5, label = 'pruning')
plt.ylabel("# of event")
plt.xlabel('Num_triggered')
plt.legend(loc=0)





###################################################################################
####################### Part in which the various rates are computed ##############
###################################################################################


# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
#######    enerbins = np.unique(ev_select[:,1])
#zenbins = 180-np.unique(A_rect[:,3])
# zenbins = np.array([94.77,95.74,97.18,98.21,99.59,101.54, 104.48, 106.6, 109.47, 113.58,120,132])
# zenbins = 180. - zenbins
# #zenbins = [94,100,105,110,120,131]
# stepbins = np.unique(ev_select[:,2])


m, v = ua.compute_meanNtrig(
    lay1.energy_bins_limits,
    lay1.zenith_bins_limits,
    lay1.ev_select
)


# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins
ua.plot_meanNtrig(
    m,
    v,
    lay1.energy_bins_limits,
    lay1.energy_bins_centers,
    lay1.zenith_bins_limits,
    lay1.zenith_bins_centers,
    layout=lay1.mask_name 
)


trig_rate = ua.compute_trig_rate(
    lay1.energy_bins_limits,
    lay1.zenith_bins_limits,
    lay1.ev_select
)


ua.plot_trigrate(
    trig_rate,
    lay1.energy_bins_limits,
    lay1.energy_bins_centers,
    lay1.zenith_bins_limits,
    lay1.zenith_bins_centers,
    layout=lay1.mask_name 
)



# meanNtrig_ener1, varNtrig_ener1 = ua.compute_meanNtrig(
#     stepbins,
#     enerbins,
#     zenbins,
#     ev_select_trihex_nopruning
# )



# trigger rate calculation over full array
# convolving with measured CR flux

# log_E_eV = np.log10(enerbins*1.e18) -0.05
# log_E_eV = np.append(log_E_eV, log_E_eV[-1]+0.1)
# enerbins2 = 10**(log_E_eV) / 1e18 # in EeV
# delta_E = enerbins2[1:] - enerbins2[:-1]


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
    ev_select_trihex_nopruning
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




