import matplotlib.pyplot as plt
import numpy as np
import os
import diff_spec
import utils_analysis as ua

'''
analysis of simulation libraries to generate plots of triggered event rates
for various layout shapes/steps

First read all simulation files
select layout shape, primary type, trigger threshold 
store in an array to be dumped in file "events_data_dir/ev_select_*.npy"

If this npy file is already present in the events_data_dir directory, the simulation output files are not re-read.

Events are split in bins of energy, zenith angle, layout step size.
The mean and variance of the triggered number of events is calculated for a given trigger threshold.
The triggered event rate is calculated for a given setup and for given bins, 
for a minimum number of triggered antennas n_trig_thres

Plots are saved in "plot_path" directory.
'''

# read files

#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/InterpolationOutputExample/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/StshpLibrary02/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/P2PdataNew/"
path = "/Users/kotera/BROQUE/Data_GRAND/Matias/StshpLibrary-HDF5-Grids/"
#path = "/Users/benoitl/Documents/GRAND/StshpLibaryHDF5-Grids/"

plot_path = '/Users/kotera/BROQUE/Plots_GRAND/'
#plot_path = '/Users/benoitl/Documents/GRAND/plots_tests'
events_data_dir = "/Users/kotera/BROQUE/Data_GRAND/Matias/event_data"
#events_data_dir = "/Users/benoitl/Documents/GRAND/event_data"

os.makedirs(plot_path, exist_ok=True)
os.makedirs(events_data_dir, exist_ok=True)

threshold = 30 # trigger threshold for individual antennas in muV
n_trig_thres = 5 # number of triggered antennas required to trigger an event


ev_select_rect_file = os.path.join(events_data_dir, 'ev_select_rect_proton.npy')
ev_select_hexhex_file = os.path.join(events_data_dir, 'ev_select_hexhex_proton.npy')
ev_select_trihex_file = os.path.join(events_data_dir, 'ev_select_trihex_proton.npy')


## Check that A_rect exits

do_make_ev_list = (
    ua.isnot_ev_select(ev_select_rect_file) *
    ua.isnot_ev_select(ev_select_hexhex_file) *
    ua.isnot_ev_select(ev_select_trihex_file)
)

if do_make_ev_list:
    ev_list = ua.make_ev_list(path)
    # select triggered antennas and events
    
    for ev in ev_list:
        if "voltage" in ev.name:
            ev.num_triggered = sum(ev.is_triggered1(threshold))
            ev.is_triggered2 = (ev.num_triggered > n_trig_thres)


    ev_select_rect = ua.make_ev_select(
        ev_list,
        'rect',
        'Proton', # 'Proton', 'Fe^56', 'Gamma'
        ev_select_rect_file
    )

    ev_select_hexhex = ua.make_ev_select(
        ev_list,
        'hexhex',
        'Proton',
        ev_select_hexhex_file
    )
    ev_select_trihex = ua.make_ev_select(
        ev_list,
        'trihex',
        'Proton',
        ev_select_trihex_file
    )

else:
    ev_select_rect = np.load(ev_select_rect_file)
    ev_select_hexhex = np.load(ev_select_hexhex_file)
    ev_select_trihex = np.load(ev_select_trihex_file)


ev_select = ev_select_rect
layout = "rect"

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
enerbins = np.unique(ev_select[:,1])
#zenbins = 180-np.unique(A_rect[:,3])
zenbins = np.array([94.77,95.74,97.18,98.21,99.59,101.54, 104.48, 106.6, 109.47, 113.58,120,132])
zenbins = 180. - zenbins
#zenbins = [94,100,105,110,120,131]
stepbins = np.unique(ev_select[:,2])

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)




# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)



sym_list = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
myc = ['0','0.20','0.4','0.6','0.8']



do_plot_rate_fixedenergy_vszenith = True
do_plot_rate_fixedzenith_vsenergy = False
do_plot_rate_fixedsteps_vsenergy = False
do_plot_Ntrig_fixedernergy_vszenith = False
do_plot_Ntrig_fixedzenith_vsenergy = False
do_plot_Ntrig_fixedsteps_vsenergy = False

do_plot_overall_rate_fixedenergy_vszenith = False
do_plot_overall_rate_fixedzenith_vsenergy = False
do_plot_overall_rate_fixedenergy_vszenith_allgeom = False
do_plot_overall_rate_fixedzenith_vsenergy_allgeom = False
do_plot_overall_rate_fixedzenith_vsenergy_allprim = False



if do_plot_Ntrig_fixedsteps_vsenergy:
# plot Ntriggered antennas vs energies for fixed steps
    ua.plot_Ntrig_fixedsteps_vsenergy(
    meanNtrig_ener,
    varNtrig_ener,
    stepbins,
    enerbins,
    zenbins,
    layout
)

# plot Ntriggered antennas vs energies for fixed zenith angles
if do_plot_Ntrig_fixedzenith_vsenergy:
    ua.plot_Ntrig_fixedzenith_vsenergy(
    meanNtrig_ener,
    varNtrig_ener,
    stepbins,
    enerbins,
    zenbins,
    layout
)


if do_plot_Ntrig_fixedernergy_vszenith:
    ua.plot_Ntrig_fixedernergy_vszenith(
    meanNtrig_ener,
    varNtrig_ener,
    stepbins,
    enerbins,
    zenbins,
    layout,
    plot_path=plot_path
)


if do_plot_rate_fixedsteps_vsenergy:
    # plot Ntriggered events vs energies for fixed step size and zenith angles
    ua.plot_rate_fixedsteps_vsenergy(
    Ntrig2_ener,
    stepbins,
    enerbins,
    zenbins,
    layout,
    plot_path=plot_path
)

if do_plot_rate_fixedzenith_vsenergy:
    ua.plot_rate_fixedzenith_vsenergy(
    Ntrig2_ener,
    stepbins,
    enerbins,
    zenbins,
    layout,
    plot_path=plot_path    
)


if do_plot_rate_fixedenergy_vszenith:
    ua.plot_rate_fixedenergy_vszenith(
    Ntrig2_ener,
    stepbins,
    enerbins,
    zenbins,
    layout,
    plot_path=plot_path    
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

area = stepbins**2 * 196 

rate = Ntrig2_ener.copy() * 0
for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            rate[istep, iener, izen] = (
                Ntrig2_ener[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area[istep] * np.cos(zen*np.pi/180) / area[istep] # !! finally we divide by the area here
            )


if do_plot_overall_rate_fixedzenith_vsenergy:          
    for izen in range(0, len(zenbins)-1):
        plt.figure(izen) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                enerbins,
                rate[istep,:,izen] * 24*3600,
                fmt=sym_list[istep],
                ls='-',
                capsize=2,
                alpha=0.7,
                label='step = %d m'%(np.int32(step))
            )
        plt.yscale('log')
        plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} m^{-2}]$')
        plt.xlabel('energy [EeV]')
        plt.title('%s, %4.2f > zenith > %4.2f deg'%(layout, thetar[izen], thetal[izen+1]))
        plt.legend(loc=1)
        #plt.ylim(1.e-1,1.e2)
        plt.ylim(1.e-11,1.e-5)
        #plt.show()
        plt.savefig(os.path.join(plot_path,'evrate_vs_energy_z%4.1f_%s_30muV.png'%(180-zenbins[izen+1], layout)))


if do_plot_overall_rate_fixedenergy_vszenith:
    for iener, ener in enumerate(enerbins):
        plt.figure(iener) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                zenbins,
                rate[istep,iener,:] * 24*3600 ,  
                fmt=sym_list[istep],
                ls='-',
                capsize=2,
                alpha=0.7,
                label='step = %d m'%(np.int32(step))
            )
        plt.yscale('log')
        plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} m^{-2}]$')
        plt.xlabel('zenith [deg]')
        plt.title('Proton, %s, E = %4.3f EeV'%(layout, ener))
        plt.legend(loc=1)
        #plt.ylim(1.e-1,1.e2)
        plt.ylim(1.e-11,1.e-5)
        plt.xlim(45,90)
        #plt.show()
        plt.savefig(os.path.join(plot_path, 'evrate_vs_zen_E%4.3f_%s_Proton.png'%(ener, layout)))





# triggered event rate (per day) over the full area, 
# convolved with the measured cosmic-ray flux
# plot all 3 geometrical layouts in one single plot

zenbins = np.array([94.77, 95.74, 97.18, 98.21, 99.59, 101.54, 104.48, 106.6, 109.47, 113.58, 120, 132])
zenbins = 180. - zenbins
enerbins = np.unique(ev_select[:,1])

ev_select = ev_select_rect

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
stepbins = np.unique(ev_select[:,2])
stepbins_rect = stepbins

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_rect = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

ev_select = ev_select_hexhex

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
stepbins = np.unique(ev_select[:,2])
stepbins_hexhex = stepbins

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_hexhex = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

ev_select = ev_select_trihex

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
stepbins = np.unique(ev_select[:,2])
stepbins_trihex = stepbins

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_trihex = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
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

area_rect = stepbins_rect**2 * 196 
area_hexhex = stepbins_hexhex**2 *3*np.sqrt(3)/2 * 91
area_trihex = stepbins_trihex**2 *np.sqrt(3)/4 * 61 * 6

rate_rect = Ntrig2_ener_rect.copy() * 0
rate_hexhex = Ntrig2_ener_hexhex.copy() * 0
rate_trihex = Ntrig2_ener_trihex.copy() * 0

for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            rate_rect[istep, iener, izen] = (
                Ntrig2_ener_rect[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_rect[istep] * np.cos(zen*np.pi/180) / area_rect[istep] # !! finally we divide by the area here
            )
            rate_hexhex[istep, iener, izen] = (
                Ntrig2_ener_hexhex[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_hexhex[istep] * np.cos(zen*np.pi/180) /area_hexhex[istep] # !! finally we divide by the area here
            )
            rate_trihex[istep, iener, izen] = (
                Ntrig2_ener_trihex[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_trihex[istep] * np.cos(zen*np.pi/180) /area_trihex[istep] # !! finally we divide by the area here
            )


if do_plot_overall_rate_fixedzenith_vsenergy_allgeom:
    for izen in range(0, len(zenbins)-1):
        plt.figure(izen) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                enerbins,
                rate_rect[istep,:,izen] * 24*3600,
                fmt="C%d"%istep+'s', 
                ms=7,
                ls='-',
                capsize=2,
                alpha=0.7,
                label='step = %d m'%(np.int32(step))
            )
            plt.errorbar(
                enerbins,
                rate_hexhex[istep,:,izen] * 24*3600,
                fmt="C%d"%istep+'h',
                ms=10,           
                ls='-',
                capsize=2,
                alpha=0.7,
                #label='step = %d m'%(np.int32(step))
            )
            plt.errorbar(
                enerbins,
                rate_trihex[istep,:,izen] * 24*3600,
                fmt="C%d"%istep+'^',
                ms=7,
                ls='-',
                capsize=2,
                alpha=0.7,
                #label='step = %d m'%(np.int32(step))
            )
        plt.yscale('log')
        plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} m^{-2}]$')
        plt.xlabel('energy [EeV]')
        plt.title('%s, %4.2f > zenith > %4.2f deg'%(layout, thetar[izen], thetal[izen+1]))
        plt.legend(loc=1)
        #plt.ylim(1.e-1,1.e2)
        plt.ylim(1.e-11,1.e-5)        
        #plt.show()
        plt.savefig(os.path.join(plot_path,'evrate_vs_energy_z%4.1f_all_30muV.png'%(180-zenbins[izen+1])))



# plots for three primaries

# assume these files have been produced first using first part of this script
ev_select_rect_file = os.path.join(events_data_dir, 'ev_select_rect_proton.npy')
ev_select_rect_fe_file = os.path.join(events_data_dir, 'ev_select_rect_fe.npy')
ev_select_rect_gamma_file = os.path.join(events_data_dir, 'ev_select_rect_gamma.npy')

ev_select_rect = np.load(ev_select_rect_file)
ev_select_rect_fe = np.load(ev_select_rect_fe_file)
ev_select_rect_gamma = np.load(ev_select_rect_gamma_file)

#####protons
ev_select = ev_select_rect

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
stepbins = np.unique(ev_select[:,2])
stepbins_rect = stepbins

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_proton = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

#####iron
ev_select = ev_select_rect_fe

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_fe = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

#####gamma
ev_select = ev_select_rect_gamma

# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 

meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select
)

# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener_gamma = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select
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

area_rect = stepbins_rect**2 * 196 

rate_proton = Ntrig2_ener_proton.copy() * 0
rate_fe = Ntrig2_ener_fe.copy() * 0
rate_gamma = Ntrig2_ener_gamma.copy() * 0

for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            rate_proton[istep, iener, izen] = (
                Ntrig2_ener_proton[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_rect[istep] * np.cos(zen*np.pi/180)/area_rect[istep]
            )
            rate_fe[istep, iener, izen] = (
                Ntrig2_ener_fe[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_rect[istep] * np.cos(zen*np.pi/180) /area_rect[istep]
            )
            rate_gamma[istep, iener, izen] = (
                Ntrig2_ener_gamma[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area_rect[istep] * np.cos(zen*np.pi/180)/area_rect[istep]
            )

colors = ['k', 'r', 'g', 'c']

if do_plot_overall_rate_fixedzenith_vsenergy_allprim:
    for izen in range(0, len(zenbins)-1):
        plt.figure(izen) 
        plt.clf()
        plot_lines = []
        for istep, step in enumerate(stepbins):
            l1, = plt.plot(
                enerbins,
                rate_proton[istep,:,izen] * 24*3600,
                "C%d"%istep+'^',
                ms=7,
                ls='-',
                alpha=0.7,
                #label='step = %d m'%(np.int32(step))
            )
            l2, = plt.plot(
                enerbins,
                rate_fe[istep,:,izen] * 24*3600,
                "C%d"%istep+'s',
                ms=10,           
                ls='-',
                alpha=0.7,
            )
            l3, = plt.plot(
                enerbins,
                rate_gamma[istep,:,izen] * 24*3600,
                "C%d"%istep+'o',
                ms=7,
                ls='-',
                alpha=0.7,
            )
            plot_lines.append([l1, l2, l3])

        plt.yscale('log')
        plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1} m^{-2}]$')
        plt.xlabel('energy [EeV]')
        plt.title('%s, %4.2f > zenith > %4.2f deg'%(layout, thetar[izen], thetal[izen+1]))
        #plt.legend(loc=1)
        legend1 = plt.legend(plot_lines[0], ["Proton", "Fe", "Gamma"], loc=2)
        parameters = ['step = %d m'%(np.int32(step)) for step in stepbins]
        plt.legend([l[0] for l in plot_lines], parameters, loc=1)
        plt.gca().add_artist(legend1)
        #plt.legend([lines[i] for i in [0,1,2]], ["Proton", "Fe", "Gamma"], loc=2)
        #plt.ylim(1.e-1,1.e2)
        plt.ylim(1.e-10,1.e-4)
        #plt.show()
        plt.savefig(os.path.join(plot_path,'evrate_vs_energy_z%4.1f_allprim_30muV.png'%(180-zenbins[izen+1])))


# geometrical factor
dtheta = thetal - thetar
plt.figure()
plt.clf()
plt.plot(zenbins, np.cos((90-zenbins)*np.pi/180) * np.sin((90-zenbins)*np.pi/180) * dtheta)
plt.xlabel('zenith angle')
plt.ylabel("angular factor")
plt.yscale('log')
plt.savefig(os.path.join(plot_path,'angular_factor.png'))

