import numpy as np
import os
import matplotlib.pyplot as plt

'''
Definition of class to read the simulation output files
and of event selection routines 
to perform simulation output analysis in the analysis.py routine
'''

SYM_LIST = ['.','o','v','*','s','.','o','v','*','s','.','o','v','*','s']
MYC = ['0','0.20','0.4','0.6','0.8']


class Event:
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


def isnot_ev_select(ev_select_file):
    return not(os.path.isfile(ev_select_file))


def make_ev_select(ev_list, layout, primary, ev_select_file):
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


def compute_meanNtrig(stepbins, enerbins, zenbins, ev_select):

    meanNtrig_ener = []
    varNtrig_ener = []

    for istep, step in enumerate(stepbins):
        meanNtrig_step = []
        varNtrig_step = []

        for iener, ener in enumerate(enerbins):
            meanNtrig_zen = []
            varNtrig_zen = []
        
            #for izen in range(0, len(zenbins)-1):
            for izen, zen in enumerate(zenbins):
                ind = np.where((ev_select[:,1] == ener) * (ev_select[:,2] == step)
                    * (np.abs(ev_select[:,3]-(180-zen)) < 0.5) * (ev_select[:,0] > 0))
                    #* (A_rect[:,3] == 180-zen)* (A_rect[:,0] > 0))
                    #* (A_rect[:,3] >= zenbins[izen]) * (A_rect[:,3] < zenbins[izen+1]))
                #print(ind)
                if (len(ind[0]) == 0):
                    meanNtrig_zen.append(0)
                    varNtrig_zen.append(0)   
                else:
                    meanNtrig_zen.append(np.mean(ev_select[ind[0],0]))
                    varNtrig_zen.append(np.var(ev_select[ind[0],0]))

            meanNtrig_step.append(meanNtrig_zen)
            varNtrig_step.append(varNtrig_zen)

        meanNtrig_ener.append(meanNtrig_step)
        varNtrig_ener.append(varNtrig_step)

    meanNtrig_ener = np.array(meanNtrig_ener)
    varNtrig_ener = np.array(varNtrig_ener)
    return meanNtrig_ener, varNtrig_ener


def compute_trig_rate(stepbins, enerbins, zenbins, ev_select):

    Ntrig2_ener = []

    for istep, step in enumerate(stepbins):
        Ntrig2_step = []

        for iener, ener in enumerate(enerbins):
            Ntrig2_zen = []
        
            #for izen in range(0, len(zenbins)-1):
            for izen, zen in enumerate(zenbins):
                ind = np.where((ev_select[:,1] == ener) * (ev_select[:,2] == step) 
                    *(np.abs(ev_select[:,3]-(180-zen)) < 0.5))
                    #* (A_rect[:,3] == 180-zen))
                    #* (A_rect[:,3] >= zenbins[izen]) * (A_rect[:,3] < zenbins[izen+1]))
                if len(ind[0])==0:
                    Ntrig2_zen.append(0)
                else:
                    Ntrig2_zen.append(sum(ev_select[ind[0],4])/np.size(ev_select[ind[0],4]))

            Ntrig2_step.append(Ntrig2_zen)

        Ntrig2_ener.append(Ntrig2_step)

    Ntrig2_ener = np.array(Ntrig2_ener)
    return Ntrig2_ener


def plot_Ntrig_fixedsteps_vsenergy(
    meanNtrig_ener,
    varNtrig_ener,
    stepbins,
    enerbins,
    zenbins
):
    """
    plot Ntriggered antennas vs energies for fixed steps
    """
    for istep, step in enumerate(stepbins):
        plt.figure(istep) 
        plt.clf()
        for izen in range(0, len(zenbins)-1):
            plt.errorbar(
                enerbins,
                meanNtrig_ener[istep,:,izen],
                yerr=np.sqrt(varNtrig_ener[istep,:,izen]), 
                fmt=SYM_LIST[izen],
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


def plot_Ntrig_fixedzenith_vsenergy(
    meanNtrig_ener,
    varNtrig_ener,
    stepbins,
    enerbins,
    zenbins
):
    """
    plot Ntriggered antennas vs energies for fixed zenith angles
    """
    for izen in range(0, len(zenbins)-1):
        plt.figure(izen+4) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                enerbins,
                meanNtrig_ener[istep,:,izen],
                yerr=np.sqrt(varNtrig_ener[istep,:,izen]), 
                fmt=SYM_LIST[istep],
                capsize=2,
                alpha=0.7,
                label='step = %d m'%(np.int32(step)))
            #plt.errorbar(enerbins, Ntrig2_ener[istep,:,izen],  
            #   fmt=sym_list[istep], capsize=2, alpha=0.7)
        plt.yscale('log')
        plt.ylabel('N triggered antennas')
        plt.xlabel('energy [EeV]')
        plt.title('hex, %4.0f > zenith >%4.0f deg'%(180-zenbins[izen], 180-zenbins[izen+1]))
        plt.legend(loc=4)
        plt.show()


def plot_Ntrig_fixedernergy_vszenith(
    meanNtrig_ener,
    varNtrig_ener,
    stepbins,
    enerbins,
    zenbins,
    plot_path='./'
):
    """
    # plot Ntriggered antennas vs zenith angles for fixed energies 
    """
    for iener, ener in enumerate(enerbins):
        plt.figure(iener) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                zenbins,
                meanNtrig_ener[istep,iener,:],
                yerr=np.sqrt(varNtrig_ener[istep,iener,:]), 
                fmt=SYM_LIST[istep],
                capsize=2,
                alpha=0.7,
                label='step = %d m'%(np.int32(step))
            )
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


def plot_rate_fixedsteps_vsenergy(
    Ntrig2_ener,
    stepbins,
    enerbins,
    zenbins,
    plot_path='./'
):
    for istep, step in enumerate(stepbins):
        plt.figure(istep) 
        plt.clf()
        for izen in range(0, len(zenbins)-1):
            plt.errorbar(
                enerbins,
                Ntrig2_ener[istep,:,izen], 
                fmt=SYM_LIST[izen],
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


def plot_rate_fixedzenith_vsenergy(
    Ntrig2_ener,
    stepbins,
    enerbins,
    zenbins,
    plot_path="./"   
):
    for izen in range(0, len(zenbins)-1):
        plt.figure(izen+4) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                enerbins,
                Ntrig2_ener[istep,:,izen],  
                fmt=SYM_LIST[istep],
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


def plot_rate_fixedenergy_vszenith(
    Ntrig2_ener,
    stepbins,
    enerbins,
    zenbins,
    plot_path="./"  
):
# plot Ntriggered events vs zenith angles for fixed energies 
    for iener, ener in enumerate(enerbins):
        plt.figure(iener) 
        plt.clf()
        for istep, step in enumerate(stepbins):
            plt.errorbar(
                zenbins,
                Ntrig2_ener[istep,iener,:],  
                fmt=SYM_LIST[istep],
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
