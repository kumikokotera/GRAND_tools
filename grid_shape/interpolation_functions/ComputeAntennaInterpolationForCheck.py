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
import numpy as np
import matplotlib.pyplot as plt

import hdf5fileinout as hdf5io
#from matplotlib.pyplot import cm
import DatabaseFunctions as mydatabase  #my database handling library
import AiresInfoFunctions as AiresInfo

from StarshapeInterpolation import do_interpolation_hdf5


parser = argparse.ArgumentParser(description='A script to get the CPU time in a library of Simulations')
parser.add_argument('DatabaseFile', #name of the parameter
                    metavar="filename", #name of the parameter value in the help
                    help='The Database of the library .db file') # help message for this parameter
results = parser.parse_args()
dbfile=results.DatabaseFile

dbfile="/home/mjtueros/GRAND/GP300/HDF5StshpLibrary/StshpXmaxLibraryInExa24.01.sql3.db"
dbfile="/home/mjtueros/AiresRepository/DiscreteLibraryTest/ConicalStshpTestLibrary.sql3.db"
dbfile="/home/mjtueros/AiresRepository/DiscreteLibrary/StshpLibrary-01.sql3.db"
#directory where the files from the library are located
Directory = "/home/mjtueros/GRAND/GP300/HDF5StshpLibrary/Outbox"
Directory = "/home/mjtueros/AiresRepository/DiscreteLibraryTest/ConicalStshpOutbox"
Directory = "/home/mjtueros/AiresRepository/DiscreteLibrary/StshpOutbox"
#what to use in the interpolation (efield, voltage, filteredvoltage, all)
# all, efield, voltage, filtered voltage
usetrace='all'
DISPLAY=True

##logging.debug('This is a debug message')
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

while(DatabaseRecord!=None): #500 events in 30min, withouth tresholding, 700 en 47 min
    #while(DatabaseRecord!=None):
    DatabaseStatus = mydatabase.GetStatusFromRecord(DatabaseRecord) #i do it with a function call becouse if we change the database structure we dont have to change this
    #Directory = mydatabase.GetDirectoryFromRecord(DatabaseRecord)
    JobName = mydatabase.GetNameFromRecord(DatabaseRecord)
    JobDirectory = str(Directory)+"/"+str(JobName)
    Tries = mydatabase.GetTriesFromRecord(DatabaseRecord)
    #.
    #logging.debug("Reading Job " + JobName + " which was in " + DatabaseStatus + " status ("+str(Tries)+") at " + Directory)
    #.
    if(DatabaseStatus == "RunOK" and "Gamma" in JobName): # and "_3.9" in JobName and "Proton" in JobName):
        try:
            InputFilename=str(Directory)+"/"+str(JobName)+"/"+str(JobName)+".hdf5"
            print("InputFileName:"+InputFilename)
            #
            OutputFilename = InputFilename + '.Interpolated.'+str(usetrace)+'.hdf5'
            #
            if(os.path.isfile(OutputFilename)):
             print("already computed")
             DatabaseRecord=CurDataBase.fetchone()
             countok += 1
             continue
            #.
            CurrentRunInfo=hdf5io.GetRunInfo(InputFilename)
            CurrentEventName=hdf5io.GetEventName(CurrentRunInfo,0) #using the first event of each file (there is only one for now)
            CurrentAntennaInfo=hdf5io.GetAntennaInfo(InputFilename,CurrentEventName)
            #.
            #one way of putting the antenna information as the numpy array this script was designed to use:
            antennamin=160
            antennamax=176 # WE ARE GETTING THE RANDOM antennas!
            AntID=CurrentAntennaInfo['ID'].data[antennamin:antennamax]
            xpoints=CurrentAntennaInfo['X'].data[antennamin:antennamax]
            ypoints=CurrentAntennaInfo['Y'].data[antennamin:antennamax]
            zpoints=CurrentAntennaInfo['Z'].data[antennamin:antennamax]
            t0points=CurrentAntennaInfo['T0'].data[antennamin:antennamax]
            #
            NewPos=np.stack((xpoints,ypoints,zpoints,t0points), axis=1)
            #
            #
            do_interpolation_hdf5(NewPos, InputFilename, OutputFilename, antennamin=0, antennamax=159, EventNumber=0, DISPLAY=DISPLAY, usetrace=usetrace)
            #
            countok += 1
            #
            #now, lets open the Antenna
            print("Event #{} done".format(countok))
            #.
            #.
        except FileNotFoundError:
          logging.error("ant_interpol_chk_db:file not found or invalid:"+TaskName)
          counterr += 1
          #.
    else:
     print(JobName + " is in " + DatabaseStatus + " Status, Skipping")
    #this is the last order of the while, that will fetch the next record of the database
    DatabaseRecord=CurDataBase.fetchone()
#.
#.

