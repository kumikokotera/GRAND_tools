'''
Performs an interpolation of the peak-to-peak electric field at 6 antenna positions for each shower in database.
Check interpolation efficiency.

run ant_trig_db [argv1]

argv1: str
    path+name of shower database
'''
import sys
import os
import logging   #for...you guessed it...logging
import sqlite3   #for the database
import argparse  #for command line parsing
import glob      #for listing files in directories
import importlib #to be able to import modules dynamically (will need it to switch from cluster to local run configurations)
import time      #for the sleep
import datetime  #for the now()
#sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/Scripts") #so that it knows where to find things
import numpy as np
import matplotlib.pyplot as plt
import hdf5fileinout as hdf5io
#from matplotlib.pyplot import cm
import DatabaseFunctions as mydatabase  #my database handling library
import AiresInfoFunctions as AiresInfo
#import interpol_func as intf
import interpol_func_hdf5 as intf
import grids

parser = argparse.ArgumentParser(description='A script to get the CPU time in a library of Simulations')
parser.add_argument('DatabaseFile', #name of the parameter
                    metavar="filename", #name of the parameter value in the help
                    help='The Database of the library .db file') # help message for this parameter
results = parser.parse_args()
dbfile=results.DatabaseFile

dbfile="/home/mjtueros/GRAND/GP300/HDF5StshpLibrary/StshpXmaxLibraryInExa24.01.sql3.db"
#directory where the files from the library are located
Directory = "/home/mjtueros/GRAND/GP300/HDF5StshpLibrary/Outbox"
#what to use in the interpolation (efield, voltage, filteredvoltage)
usetrace='efield'
#threshold abouve wich the interpolation is computed
threshold=0;#26 #8.66 for 15uV , 26 for 45uV
trigger=30
display=False
#logging.debug('This is a debug message')
#logging.info('This is an info message')
#logging.warning('This is a warning message')
#logging.error('This is an error message')
#logging.critical('This is a critical message')
logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.WARNING)
#logging.disable(logging.DEBUG)
logging.getLogger('matplotlib.font_manager').disabled = True


logging.debug("Starting Xmax analysis with %s " % dbfile)

DataBase=mydatabase.ConnectToDataBase(dbfile)
#this is to show the current status of the database
mydatabase.GetDatabaseStatus(DataBase)

#This is how you search on the database, here im selecting everything (To Do: functions to search the database)
#This is to get a cursor on the database. You can think of the cursor as a working environment. You can have many cursors.
CurDataBase = DataBase.cursor()
CurDataBase.execute("SELECT * FROM showers")

#f = open(dbfile+".Xmax", 'w')
#f.write("Id, Energy, Zenith, Azimuth, Primary, SlantXmax, XmaxAltitude, Xxmax, Yxmax,Zxmax\n")

DatabaseRecord = CurDataBase.fetchone()
counterr = 0
countok = 0
InterpErrAll = np.zeros((16,1100))
P2pAll = np.zeros((16,1100))
P2pAllx = np.zeros((16,1100))
P2pAlly = np.zeros((16,1100))
P2pAllz = np.zeros((16,1100))
errortype = np.zeros((16,1100))
errortypex = np.zeros((16,1100))
errortypey = np.zeros((16,1100))
errortypez = np.zeros((16,1100))
distance = np.zeros(1100)
energy= np.zeros(1100)
zenith= np.zeros(1100)

while(DatabaseRecord!=None and countok < 1000): #500 events in 30min, withouth tresholding, 700 en 47 min
    #while(DatabaseRecord!=None):
    DatabaseStatus = mydatabase.GetStatusFromRecord(DatabaseRecord) #i do it with a function call becouse if we change the database structure we dont have to change this
    #Directory = mydatabase.GetDirectoryFromRecord(DatabaseRecord)
    JobName = mydatabase.GetNameFromRecord(DatabaseRecord)
    JobDirectory = str(Directory)+"/"+str(JobName)
    Tries = mydatabase.GetTriesFromRecord(DatabaseRecord)
    #.
    logging.debug("Reading Job " + JobName + " which was in " + DatabaseStatus + " status ("+str(Tries)+") at " + Directory)
    #.
    if(DatabaseStatus == "RunOK"): #and JobName=="Stshp_XmaxLibrary_0.1995_80.40_180_Gamma_01"):
        try:
            TaskName = mydatabase.GetTasknameFromRecord(DatabaseRecord)
            Id = mydatabase.GetIdFromRecord(DatabaseRecord)
            Energy = mydatabase.GetEnergyFromRecord(DatabaseRecord)
            Zenith = mydatabase.GetZenithFromRecord(DatabaseRecord)
            Azimuth = mydatabase.GetAzimuthFromRecord(DatabaseRecord)
            Primary = mydatabase.GetPrimaryFromRecord(DatabaseRecord)
            Xmax = mydatabase.GetXmaxFromRecord(DatabaseRecord)
            #.
            InputFilename=str(Directory)+"/"+str(JobName)+"/"+str(JobName)+".hdf5"
            #.
            CurrentRunInfo=hdf5io.GetRunInfo(InputFilename)
            CurrentEventName=hdf5io.GetEventName(CurrentRunInfo,0) #using the first event of each file (there is only one for now)
            CurrentAntennaInfo=hdf5io.GetAntennaInfo4(InputFilename,CurrentEventName)
            #.
            #AntNum, AntPos, AntID = intf.get_antenna_pos_zhaires(JobDirectory+'/antpos.dat')
            #one way of putting the antenna information as the numpy array this script was designed to use:
            antennamin=160
            antennamax=176 #NOTE THAT WE ARE GETTING only THE RANDOM!
            AntID=CurrentAntennaInfo['ID'].data[antennamin:antennamax]
            #.
            #this gets the p2p values in all chanels, for all simulated antennas.
            #.
            p2p_total_sim=CurrentAntennaInfo['P2P_efield'].data[antennamin:antennamax]
            p2p_x_sim=CurrentAntennaInfo['P2Px_efield'].data[antennamin:antennamax]
            p2p_y_sim=CurrentAntennaInfo['P2Py_efield'].data[antennamin:antennamax]
            p2p_z_sim=CurrentAntennaInfo['P2Pz_efield'].data[antennamin:antennamax]
            t0_sim=CurrentAntennaInfo['HilbertPeakTime'].data[antennamin:antennamax]
            #compute the footprint size:
            p2p_all=CurrentAntennaInfo['P2P_efield'].data[0:antennamax]
            xpoints=CurrentAntennaInfo['X'].data[0:antennamax]#but for the positions we get all, cos we need all
            ypoints=CurrentAntennaInfo['Y'].data[0:antennamax]
            zpoints=CurrentAntennaInfo['Z'].data[0:antennamax]
            AntPos=np.stack((xpoints,ypoints,zpoints), axis=0)
            ind = np.where(p2p_total_sim > trigger) #now i remove the cases where the signal is less than trigger
            myAntPos= AntPos[:,ind]
            modulus = np.sqrt(myAntPos[0,:]**2+myAntPos[1,:]**2)
            if(modulus.size!=0):
              maxdistance=max(modulus[0])
            else:
              maxdistance=0
            #.
            energy[countok]=Energy
            zenith[countok]=Zenith
            distance[countok]=maxdistance
            #.
            #now, lets open the Interpolated file
            antennamin=0
            antennamax=16
            InputFilename=str(Directory)+"/"+str(JobName)+"/"+str(JobName)+".hdf5.Interpolated.efield.hdf5"
            CurrentAntennaInfo=hdf5io.GetAntennaInfo4(InputFilename,CurrentEventName)
            p2p_total_int=CurrentAntennaInfo['P2P_efield'].data[antennamin:antennamax]
            p2p_x_int=CurrentAntennaInfo['P2Px_efield'].data[antennamin:antennamax]
            p2p_y_int=CurrentAntennaInfo['P2Py_efield'].data[antennamin:antennamax]
            p2p_z_int=CurrentAntennaInfo['P2Pz_efield'].data[antennamin:antennamax]
            t0_int=CurrentAntennaInfo['HilbertPeakTime'].data[antennamin:antennamax]
            InterpErr=t0_int-t0_sim
            #
            InterpErrAll[:,countok] = InterpErr
            #.
            #and this is the p2p value of all interpolated antennas
            P2pAll[:,countok] = p2p_total_sim
            #.
            P2pAllx[:,countok] = p2p_x_int
            P2pAlly[:,countok] = p2p_y_int
            P2pAllz[:,countok] = p2p_z_int
            #.
            #now, false positives is 2, false negatives is -2, 1  is correct and triggered, -1 is neither
            a=np.arange(0,len(p2p_total_sim))
            errortype[:,countok]=[2 if  (p2p_total_int[i] >= trigger and p2p_total_sim[i]<trigger) else -2 if (p2p_total_int[i] < trigger and p2p_total_sim[i]>=trigger) else 1 if (p2p_total_int[i] >= trigger and p2p_total_sim[i]>=trigger) else -1  for i in a]
            errortypex[:,countok]=[2 if  (p2p_x_int[i] >= trigger and p2p_x_sim[i]<trigger) else -2 if (p2p_x_int[i] < trigger and p2p_x_sim[i]>=trigger) else 1 if (p2p_x_int[i] >= trigger and p2p_x_sim[i]>=trigger) else -1 for i in a]
            errortypey[:,countok]=[2 if  (p2p_y_int[i] >= trigger and p2p_y_sim[i]<trigger) else -2 if (p2p_y_int[i] < trigger and p2p_y_sim[i]>=trigger) else 1 if (p2p_y_int[i] >= trigger and p2p_y_sim[i]>=trigger) else -1 for i in a]
            errortypez[:,countok]=[2 if  (p2p_z_int[i] >= trigger and p2p_z_sim[i]<trigger) else -2 if (p2p_z_int[i] < trigger and p2p_z_sim[i]>=trigger) else 1 if (p2p_z_int[i] >= trigger and p2p_z_sim[i]>=trigger) else -1 for i in a]
            #.
            countok += 1
            #
            print("Event #{} done".format(countok))
            #.
        except FileNotFoundError:
          logging.error("ant_interpol_chk_db:file not found or invalid:"+TaskName)
          counterr += 1
          #.
    #this is the last order of the while, that will fetch the next record of the database
    DatabaseRecord=CurDataBase.fetchone()
#.
#.
logging.debug('ant_interpol_chk_db: Plotting...')

ind = np.where(P2pAll != 0) #now i remove the cases where the signal is 0
myInterpErrAll = InterpErrAll[ind]
myP2pAll= P2pAll[ind]

errortype=errortype[:,0:countok] #here i remove the unused part of the arrays
errortypex=errortypex[:,0:countok]
errortypey=errortypey[:,0:countok]
errortypez=errortypez[:,0:countok]

errortype=errortype.flatten()
errortypex=errortypex.flatten()
errortypey=errortypey.flatten()
errortypez=errortypez.flatten()

##############Plot histogram of relative errors, for all components####################################
fig1 = plt.figure(1,figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
mybins = [-2.5,-1.5,-0.5,0.5,1.5,2.5]

ax1=fig1.add_subplot(221)
ax1.set_xlabel('type:-2 false negative -1 negative 1 positive 2 false positive')
ax1.set_ylabel('N')
name = 'clasification errors ' + str(usetrace) + " threshold " + str(threshold) + " trigger " + str(trigger)
plt.title(name)
plt.yscale('log')
plt.hist(errortype, bins=mybins,alpha=0.8,label="Total",density=True)

ax2=fig1.add_subplot(222)
ax2.set_xlabel('type:-2 false negative -1 negative 1 positive 2 false positive')
ax2.set_ylabel('N')
name = 'clasification errors x' + str(usetrace) + " threshold " + str(threshold) + " trigger " + str(trigger)
plt.title(name)
plt.yscale('log')
plt.hist(errortypex, bins=mybins,alpha=0.8,label="Total",density=True)

ax3=fig1.add_subplot(223)
ax3.set_xlabel('type:-2 false negative -1 negative 1 positive 2 false positive')
ax3.set_ylabel('N')
name = 'clasification errors y' + str(usetrace) + " threshold " + str(threshold) + " trigger " + str(trigger)
plt.title(name)
plt.yscale('log')
plt.hist(errortypey, bins=mybins,alpha=0.8,label="Total",density=True)

ax4=fig1.add_subplot(224)
ax4.set_xlabel('type:-2 false negative -1 negative 1 positive 2 false positive')
ax4.set_ylabel('N')
name = 'clasification errors z ' + str(usetrace) + " threshold " + str(threshold) + " trigger " + str(trigger)
plt.title(name)
plt.yscale('log')
plt.hist(errortypez, bins=mybins,alpha=0.8,label="Total",density=True)

plt.show()
##############Plot histogram of relative errors, for all components####################################
fig2 = plt.figure(2,figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
mybins = np.linspace(-495,495,199)

ax1=fig2.add_subplot(111)
ax1.set_xlabel('$t_{int}-t_{sim}|$')
ax1.set_ylabel('N')
name = 'overall errors ' + str(usetrace) + " threshold " + str(threshold)
plt.title(name)
plt.hist(myInterpErrAll, bins=mybins,alpha=0.8,label="Total")
plt.tight_layout()

####################Plot 2d histogram, relative errors vs signal, all components##################################3

fig3 = plt.figure(3,figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
ax1=fig3.add_subplot(111)
name = ' Total ' + str(usetrace) + " threshold " + str(threshold)
plt.title(name)
ax1.set_ylabel('$t_{int}-t_{sim}$')
ax1.set_xlabel('$log_{10} E_{sim}  [\mu V/m]$')
plt.hist2d(np.log10(myP2pAll),myInterpErrAll,bins=[100,200],range=[[-2, 3], [-50, 50]])

##########################Plot scatter 2d errors vs signal all components###########################3

fig4 = plt.figure(4,figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
ax1=fig4.add_subplot(111)
name = 'scatter overall errors'
plt.title(name)
ax1.set_ylabel('$t_{int}-t_{sim} (ns)$')
ax1.set_xlabel('$log_{10} E_{sim}  [\mu V/m]$')

plt.scatter(np.log10(myP2pAll),myInterpErrAll, s=1)
plt.xlim(-2,3)
plt.ylim(-495, 495)

plt.tight_layout()

fig5 = plt.figure(5,figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
ax1=fig5.add_subplot(211)
name = 'size vs energy'
plt.title(name)
ax1.set_ylabel('size')
ax1.set_xlabel('energy$')
plt.scatter(energy,distance, s=1)
#plt.xlim(-2,3)
#plt.ylim(-495, 495)

ax2=fig5.add_subplot(212)
name = 'size vs zenith'
plt.title(name)
ax1.set_ylabel('size')
ax1.set_xlabel('zenith$')
plt.scatter(180-zenith,distance, s=1)
plt.xlim(40,90)

plt.tight_layout()



plt.show()
