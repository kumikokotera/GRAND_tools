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
for a minimum number of triggered antennas Ntrig_thres

Plots are saved in "plot_path" directory.
'''

# read files

#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/InterpolationOutputExample/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/StshpLibrary02/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/P2PdataNew/"
#path = "/Users/kotera/BROQUE/Data_GRAND/Matias/StshpLibaryHDF5-Grids/"
path = "/Users/benoitl/Documents/GRAND/StshpLibaryHDF5-Grids/"

#plot_path = '/Users/kotera/BROQUE/Plots_GRAND/'
plot_path = '/Users/benoitl/Documents/GRAND/plots_tests'
events_data_dir = "/Users/benoitl/Documents/GRAND/event_data"

os.makedirs(plot_path, exist_ok=True)
os.makedirs(events_data_dir, exist_ok=True)

threshold = 30 # trigger threshold for individual antennas in muV
Ntrig_thres = 10 # number of triggered antennas required to trigger an event

ev_select_rect_file = os.path.join(events_data_dir, 'ev_select_rect.npy')
ev_select_hexhex_file = os.path.join(events_data_dir, 'ev_select_hexhex.npy')
ev_select_trihex_file = os.path.join(events_data_dir, 'ev_select_trihex.npy')

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
            ev.is_triggered2 = (ev.num_triggered > Ntrig_thres)


    ev_select_rect = ua.make_ev_select(
        ev_list,
        'rect',
        'Proton',
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




# calculate mean and variance of triggered antenna numbers in each zenith angle and energy bins 
enerbins = np.unique(ev_select_rect[:,1])
#zenbins = 180-np.unique(A_rect[:,3])
zenbins = np.array([94.77,95.74,97.18,98.21,99.59,101.54, 104.48, 106.6, 109.47, 113.58,120,132])
zenbins = 180. - zenbins
#zenbins = [94,100,105,110,120,131]
stepbins = np.unique(ev_select_rect[:,2])



meanNtrig_ener, varNtrig_ener = ua.compute_meanNtrig(
    stepbins,
    enerbins,
    zenbins,
    ev_select_rect
)




# calculate mean and variance of triggered event numbers in each zenith angle and energy bins 

Ntrig2_ener = ua.compute_trig_rate(
    stepbins,
    enerbins,
    zenbins,
    ev_select_rect
)


# plot Ntriggered antennas vs energies for fixed steps
sym_list = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
myc = ['0','0.20','0.4','0.6','0.8']

for istep, step in enumerate(stepbins):
    plt.figure(istep) 
    plt.clf()
    for izen in range(0, len(zenbins)-1):
        plt.errorbar(
            enerbins,
            meanNtrig_ener[istep,:,izen],
            yerr=np.sqrt(varNtrig_ener[istep,:,izen]), 
            fmt=sym_list[izen],
            capsize=2,
            alpha=0.7,
            label='%4.0f > zen >%4.0f deg'%(180-zenbins[izen],180-zenbins[izen+1])
        )
        #plt.errorbar(enerbins, Ntrig2_ener[istep,:,izen], 
        #    fmt=sym_list[izen], capsize=2)
    plt.yscale('log')
    plt.ylabel('N triggered antennas')
    plt.xlabel('energy [EeV]')
    plt.title('trihex, step = %d m'%(np.int32(step)))
    plt.legend(loc=4)
    plt.show()
 
# plot Ntriggered antennas vs energies for fixed zenith angles
for izen in range(0, len(zenbins)-1):
    plt.figure(izen+4) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        plt.errorbar(enerbins, meanNtrig_ener[istep,:,izen], yerr=np.sqrt(varNtrig_ener[istep,:,izen]), 
            fmt=sym_list[istep], capsize=2, alpha=0.7, label='step = %d m'%(np.int32(step)))
        #plt.errorbar(enerbins, Ntrig2_ener[istep,:,izen],  
         #   fmt=sym_list[istep], capsize=2, alpha=0.7)
    plt.yscale('log')
    plt.ylabel('N triggered antennas')
    plt.xlabel('energy [EeV]')
    plt.title('hex, %4.0f > zenith >%4.0f deg'%(180-zenbins[izen], 180-zenbins[izen+1]))
    plt.legend(loc=4)
    plt.show()

# plot Ntriggered antennas vs zenith angles for fixed energies 
for iener, ener in enumerate(enerbins):
    plt.figure(iener) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        plt.errorbar(zenbins, meanNtrig_ener[istep,iener,:], yerr=np.sqrt(varNtrig_ener[istep,iener,:]), 
            fmt=sym_list[istep], capsize=2, alpha=0.7, label='step = %d m'%(np.int32(step)))
        #plt.errorbar(enerbins, Ntrig2_ener[istep,:,izen], 
        #    fmt=sym_list[izen], capsize=2)
    plt.yscale('log')
    plt.ylabel('N triggered antennas')
    plt.xlabel('zenith [deg]')
    plt.title('Proton, rect, E = %4.3f EeV'%(ener))
    plt.legend(loc=2)
    plt.ylim(1,225)
    plt.xlim(45,90)
    plt.show()
    plt.savefig(os.path.join(plot_path, 'Ntrig_vs_zen_E%4.3f_rect_Proton.png'%(ener)))



# plot Ntriggered events vs energies for fixed step size and zenith angles

for istep, step in enumerate(stepbins):
    plt.figure(istep) 
    plt.clf()
    for izen in range(0, len(zenbins)-1):
        plt.errorbar(
            enerbins,
            Ntrig2_ener[istep,:,izen], 
            fmt=sym_list[izen],
            ls='-',
            capsize=2,
            alpha=0.7,
            label='%4.0f > zen >%4.0f deg'%(180-zenbins[izen], 180-zenbins[izen+1])
        )
    plt.yscale('log')
    plt.ylabel('Triggered event rate')
    plt.xlabel('energy [EeV]')
    plt.title('hex, step = %d m'%(np.int32(step)))
    plt.legend(loc=4)
    plt.show()
    plt.savefig(os.path.join(plot_path,'trigevrate_vs_energy_step%d_rect_30muV.png'%(np.int32(step))))

for izen in range(0, len(zenbins)-1):
    plt.figure(izen+4) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        plt.errorbar(
            enerbins,
            Ntrig2_ener[istep,:,izen],  
            fmt=sym_list[istep],
            ls='-',
            capsize=2,
            alpha=0.7,
            label='step = %d m'%(np.int32(step))
        )
    plt.yscale('log')
    plt.ylabel('Triggered event rate')
    plt.xlabel('energy [EeV]')
    plt.title('hex, %4.0f > zenith >%4.0f deg'%(180-zenbins[izen], 180-zenbins[izen+1]))
    plt.legend(loc=4)
    plt.show()
    plt.savefig(os.path.join(plot_path,'trigevrate_vs_energy_z%4.1f_rect_30muV.png'%(180-zenbins[izen+1])))

# plot Ntriggered events vs zenith angles for fixed energies 
for iener, ener in enumerate(enerbins):
    plt.figure(iener) 
    plt.clf()
    for istep, step in enumerate(stepbins):
        plt.errorbar(
            zenbins,
            Ntrig2_ener[istep,iener,:],  
            fmt=sym_list[istep],
            ls='-',
            capsize=2,
            alpha=0.7,
            label='step = %d m'%(np.int32(step))
        )
    plt.yscale('log')
    plt.ylabel('Triggered event rate')
    plt.xlabel('zenith [deg]')
    plt.title('Proton, rect, E = %4.3f EeV'%(ener))
    plt.legend(loc=4)
    plt.ylim(1.e-2,1.1)
    plt.xlim(45,90)
    plt.show()
    plt.savefig(os.path.join(plot_path, 'trigevrate_vs_zen_E%4.3f_rect_Proton_10N.png'%(ener)))





delta_E = enerbins[1:] - enerbins[:-1]
delta_E = np.insert(delta_E, 0, delta_E[0])

delta_omega = - (zenbins[1:] - zenbins[:-1])

delta_omega = np.insert(delta_omega, 0, delta_omega[0]) 
delta_omega = 2*np.pi * delta_omega *np.pi/180 * np.sin(np.pi/2 - zenbins*np.pi/180)

area = stepbins**2 * 200

rate = Ntrig2_ener.copy() * 0
for iener, ener in enumerate(enerbins):
    for istep, step in enumerate(stepbins):
        for izen, zen in enumerate(zenbins):
            rate[istep, iener, izen] = (
                Ntrig2_ener[istep,iener,izen] * 
                diff_spec.tale_diff_flux(ener*1e18) * 
                delta_E[iener] *1e18 * delta_omega[izen] *
                area[istep]
            )


          

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
    plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1}]$')
    plt.xlabel('energy [EeV]')
    plt.title('rect, %4.0f > zenith > %4.0f deg'%(zenbins[izen], zenbins[izen+1]))
    plt.legend(loc=4)
    plt.ylim(1.e-1,1.e2)
    plt.show()
    plt.savefig(os.path.join(plot_path,'evrate_vs_energy_z%4.1f_rect_30muV.png'%(180-zenbins[izen+1])))


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
    plt.ylabel('triggered event rate over array '+'$\\nu_{ev}\, [day^{-1}]$')
    plt.xlabel('zenith [deg]')
    plt.title('Proton, rect, E = %4.3f EeV'%(ener))
    plt.legend(loc=4)
    plt.ylim(1.e-1,1.e2)
    plt.xlim(45,90)
    plt.show()
    plt.savefig(os.path.join(plot_path, 'evrate_vs_zen_E%4.3f_rect_Proton.png'%(ener)))



















