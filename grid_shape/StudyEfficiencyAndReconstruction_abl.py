import sys
import os
import logging # for...you guessed it...logging
import sqlite3 # for the database
import argparse # for command line parsing
import glob  # for listing files in directories
import importlib # to be able to import modules dynamically (will need it to switch from cluster to local run configurations)
import time  # for the sleep
import datetime # for the now()
import numpy as np
# import scipy
import pandas as pd
import matplotlib.pyplot as plt
# ZHAIRESRUNNER = os.environ["ZHAIRESRUNNER"]
# sys.path.append(ZHAIRESRUNNER + "/TheLibraryRunner/Scripts") # so that it knows where to find things
# import DatabaseFunctions as mydatabase # my database handling library
# ZHAIRESPYTHON = os.environ["ZHAIRESPYTHON"]
# sys.path.append(ZHAIRESPYTHON)
# import AiresInfoFunctions as AiresInfo
import hdf5fileinout as hdf5io


# What i intend to do is to open all interpolated/reconstructed events, and compute:

# Detection Efficiency
# Reconstruction Efficiency

# Core Positions are Expresed with respect to FTP. The array center on the FTP is
ArrayCenter = (3400, 4800)
Xc = ArrayCenter[0]
Yc = ArrayCenter[1]
siteradius = 5250  # the circle encircling antennas of the array ()
TriggerThreshold = 45
TriggerN = 5

DetectorArea = np.pi*siteradius*siteradius
print("DetectorArea[km2]",DetectorArea/1E6)
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Creating List of files
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Library directory
# This library was run with these conditions:
# if(Id>= start and Id < = end and DatabaseStatus  == "RunOK" and Energy>1 and 180-Zenith > 50 and Primary  == "Proton"):
LibraryBaseDirectory = "/Users/ab212678/Documents/GRAND/data/OutboxAnalysis"

# since i dont have a database for these files, i make a list of all existing files, and will load everything to memmory, like in the interpolation paper
# then, i will have to see how i treat the .txt files, or if i make the effort to generate hdf5 for the lowsignal and notrigger files too.


def read_files_and_save_csv():

    EXT = "*.hdf5"
    all_hdf5_files = [
        file for path, subdir, files in os.walk(LibraryBaseDirectory) for file in glob.glob(os.path.join(path, EXT))
    ][::5]

    EXT = "*.notrigger.txt"
    all_notrigger_files = [
        file for path, subdir, files in os.walk(LibraryBaseDirectory) for file in glob.glob(os.path.join(path, EXT))
    ][::5]

    EXT = "*.toolow.txt"
    all_toolow_files = [
        file for path, subdir, files in os.walk(LibraryBaseDirectory) for file in glob.glob(os.path.join(path, EXT))
    ][::5]

    totalfiles = len(all_hdf5_files) + len(all_notrigger_files) + len(all_toolow_files)
    print(
        "hdf5 files",
        len(all_hdf5_files),
        "notrigger files",
        len(all_notrigger_files),
        "toolow files",
        len(all_toolow_files),
        "total files",
        totalfiles
    )

    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
    # Variables Definition
    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
    # EventInfo
    Zenith = np.zeros(totalfiles)
    Energy = np.zeros(totalfiles)
    Azimuth = np.zeros(totalfiles)
    CorePosition = np.zeros((totalfiles, 3))
    XmaxPosition = np.zeros((totalfiles, 3))
    NAntennas = np.zeros(totalfiles)
    NTriggered = np.zeros(totalfiles)
    # GeoRecoRef
    RecoNAntennas = np.zeros(totalfiles)
    # PlaneRecoInfo
    PlaneNAntennas = np.zeros(totalfiles)
    PlaneZenith = np.zeros(totalfiles)
    PlaneAzimuth = np.zeros(totalfiles)
    PlaneChiSignif = np.zeros(totalfiles)
    PlaneChi2 = np.zeros(totalfiles)
    # SphereRecoInfo
    SphereNAntennas = np.zeros(totalfiles)
    SphereSource = np.zeros((totalfiles, 3))
    SphereChiSignif = np.zeros(totalfiles)
    SphereChi2 = np.zeros(totalfiles)
    # ADFRecoInfo
    ADFNAntennas = np.zeros(totalfiles)
    ADFZenith = np.zeros(totalfiles)
    ADFAzimuth = np.zeros(totalfiles)
    ADFWidth = np.zeros(totalfiles)
    ADFAmplitude = np.zeros(totalfiles)
    ADFChiSignif = np.zeros(totalfiles)
    ADFChi2 = np.zeros(totalfiles)

    # SampledCores
    EventWeight = np.zeros(totalfiles)

    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
    # Start Reading files
    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 3
    loopcount = 0
    TotalEventWeight = 0
    for SimFileName in all_hdf5_files:
        if(loopcount % 100 == 0):
            print("reading", loopcount, "of", totalfiles)
        # break
        SimRunInfo = hdf5io.GetRunInfo(SimFileName)
        EventNumber = 0           # using the first event of each file (there is only one for now)
        SimEventName = hdf5io.GetEventName(SimRunInfo, EventNumber) 
        #
        # SimAntennaInfo = hdf5io.GetAntennaInfo(SimFilename,SimEventName)
        # Event Info Holds the data for this emulated event (so, it has the emulated azimuth and the emulated xmax position)
        EventInfo = hdf5io.GetEventInfo(SimFileName, SimEventName)
        Energy[loopcount] = EventInfo["Energy"][EventNumber]
        Zenith[loopcount] = EventInfo["Zenith"][EventNumber]
        Azimuth[loopcount] = EventInfo["Azimuth"][EventNumber]
        CorePosition[loopcount] = EventInfo["CorePosition"][EventNumber]
        XmaxPosition[loopcount] = EventInfo["XmaxPosition"][EventNumber]
        # AntennaP2P
        AntennaP2PInfo = hdf5io.GetAntennaP2PInfo(SimFileName, SimEventName)
        if type(AntennaP2PInfo) != type(0):
            HilbertPeakV = AntennaP2PInfo["HilbertPeakV"].data  # to get it as a numpy array
            NAntennas[loopcount] = len(HilbertPeakV)
            NTriggered[loopcount] = np.count_nonzero(HilbertPeakV > TriggerThreshold)
        else:
            # AntennaP2PCould not be found on the file. This is an error
            NAntennas[loopcount] = -3
            NTriggered[loopcount] = -3
            PlaneNAntennas[loopcount] = -3
            SphereNAntennas[loopcount] = -3 
            ADFNAntennas[loopcount] = -3
            continue
        # GeoReco
        GeoRecoRef = hdf5io.GetGeoRecoRef(SimFileName, SimEventName)
        RecoNAntennas[loopcount] = GeoRecoRef["NAntennas"]
        if((RecoNAntennas[loopcount] != NTriggered[loopcount]) and RecoNAntennas[loopcount] != -1 and TriggerThreshold == 45): # only TiggerThreshold 45 has a reconstruction
            print("RecoNAntennas is different from NTriggered...check why!",SimFileName,RecoNAntennas[loopcount],NTriggered[loopcount])
        # (becoue Ntriggered can be les than 5, and in that case reconstruction fails and thus RecoNAntennas  == -1) 
        # Plane Reconstruction
        PlaneRecoInfo = hdf5io.GetPlaneRecoInfo(SimFileName, SimEventName)
        PlaneNAntennas[loopcount] = PlaneRecoInfo["NAntennas"][EventNumber]
        PlaneZenith[loopcount] = PlaneRecoInfo["PlaneZenithRec"][EventNumber]
        PlaneAzimuth[loopcount] = PlaneRecoInfo["PlaneAzimuthRec"][EventNumber]
        PlaneChiSignif[loopcount] = PlaneRecoInfo["PlaneChiSignif"][EventNumber]
        PlaneChi2[loopcount] = PlaneRecoInfo["PlaneChi2"][EventNumber]
        # Sphere Reconstruction
        SphereRecoInfo = hdf5io.GetSphereRecoInfo(SimFileName, SimEventName)
        SphereNAntennas[loopcount] = SphereRecoInfo["NAntennas"][EventNumber]
        SphereChiSignif[loopcount] = SphereRecoInfo["SphereChiSignif"][EventNumber]
        SphereSource[loopcount] = SphereRecoInfo["SphereSource"][EventNumber]
        SphereChi2[loopcount] = SphereRecoInfo["SphereChi2"][EventNumber] 
        # ADF Reconstruction
        ADFRecoInfo = hdf5io.GetADFRecoInfo(SimFileName, SimEventName)
        ADFNAntennas[loopcount] = ADFRecoInfo["NAntennas"][EventNumber] 
        ADFZenith[loopcount] = ADFRecoInfo["ADFZenithRec"][EventNumber] 
        ADFAzimuth[loopcount] = ADFRecoInfo["ADFAzimuthRec"][EventNumber] 
        ADFWidth[loopcount] = ADFRecoInfo["WidthRec"][EventNumber] 
        ADFAmplitude[loopcount] = ADFRecoInfo["AmpRec"][EventNumber] 
        ADFChiSignif[loopcount] = ADFRecoInfo["ADFChiSignif"][EventNumber] 
        ADFChi2[loopcount] = ADFRecoInfo["ADFChi2"][EventNumber] 
        # TestedCorePositions
        ProbedCores = hdf5io.GetProbedCoresTable(SimFileName, SimEventName)
        EventWeight[loopcount] = len(ProbedCores["CoreX"])  # abl: number of time need to trigger a event
        # CoreDistances.append(np.sqrt(np.power(CorePosition[loopcount,0]-Xc,2)+np.power(CorePosition[loopcount,1]-Yc,2)) # this will include all tested core positions, with the candidate one at the end.
        #   
        TotalEventWeight += EventWeight[loopcount]
        # print(NAntennas[loopcount],NTriggered[loopcount],PlaneNAntennas[loopcount],SphereNAntennas[loopcount],ADFNAntennas[loopcount])
        loopcount += 1

    print("TotalEvents on hdf5", loopcount, "TotalWeight", TotalEventWeight)

    # temporal dataframe to get the Energy and Zenith bins
    df = pd.DataFrame(
        np.stack(
            (EventWeight, NAntennas, NTriggered, Energy, Zenith, Azimuth, PlaneNAntennas),
            axis=1), 
        columns=["EventWeight", "NAntennas", "NTriggered", "Energy", "Zenith", "Azimuth", "PlaneNAntennas"]
    )

    df.to_csv('energ_zenith_bin.csv')

    # Getting unique memebers of a Series (as the filename truncates the numbers, i have look for the true values by proximity)
    EnergyBins = df["Energy"].unique()
    EnergyBins.sort()
    ZenithBins = df["Zenith"].unique()
    ZenithBins.sort()
    ZenithBins = ZenithBins[::-1]  # reverse the order

    for i, file in enumerate(all_notrigger_files):
        if(loopcount % 100 == 0):
            print("reading", loopcount, "of", totalfiles)
        NameParts = os.path.basename(all_notrigger_files[i]).split("_")
        Primario = NameParts[1]
        Energia = float(NameParts[2])
        idx = (np.abs(EnergyBins - Energia)).argmin()
        Energia = EnergyBins[idx]
        AiresZenith = float(NameParts[3])
        idx = (np.abs(ZenithBins - (180-AiresZenith))).argmin()
        GrandZenith = ZenithBins[idx]  
        AiresAzimuth = float(NameParts[4])
        EmulatedAzimuth = float(NameParts[-1].split(".notrigger.")[0])
        Ntries = sum(1 for _ in open(all_notrigger_files[0])) # count number of lines to see how many tries there where. If this is too slow, and is constant, you can define it as a constant.
        # print(Primario,Energia,AiresZenith,GrandZenith,AiresAzimuth,EmulatedAzimuth,Ntries)
        Energy[loopcount] = Energia
        Zenith[loopcount] = GrandZenith
        Azimuth[loopcount] = EmulatedAzimuth
        EventWeight[loopcount] = Ntries
        NAntennas[loopcount] = -4
        NTriggered[loopcount] = -4
        PlaneNAntennas[loopcount] = -4
        SphereNAntennas[loopcount] = -4 
        ADFNAntennas[loopcount] = -4
        TotalEventWeight += EventWeight[loopcount] 
        loopcount += 1 

    print("TotalEvents on hdf5 + notrigger", loopcount, "TotalWeight", TotalEventWeight)

    for i, file in enumerate(all_toolow_files):
        if(loopcount % 100 == 0):
            print("reading", loopcount, "of", totalfiles)
        NameParts = os.path.basename(all_toolow_files[i]).split("_")
        Primario = NameParts[1]
        Energia = float(NameParts[2])
        idx = (np.abs(EnergyBins - Energia)).argmin()
        Energia = EnergyBins[idx]
        AiresZenith = float(NameParts[3])
        idx = (np.abs(ZenithBins - (180-AiresZenith))).argmin()
        GrandZenith = ZenithBins[idx]  
        AiresAzimuth = float(NameParts[4])
        EmulatedAzimuth = float(NameParts[-1].split(".toolow.")[0])
        # print(Primario,Energia,AiresZenith,GrandZenith,AiresAzimuth,EmulatedAzimuth,Ntries)
        Energy[loopcount] = Energia
        Zenith[loopcount] = GrandZenith
        Azimuth[loopcount] = EmulatedAzimuth
        EventWeight[loopcount] = -5
        NAntennas[loopcount] = -5
        NTriggered[loopcount] = -5
        PlaneNAntennas[loopcount] = -5
        SphereNAntennas[loopcount] = -5 
        ADFNAntennas[loopcount] = -5 
        loopcount += 1 

    print("TotalEvents on hdf5 + notrigger + toolow", loopcount, "TotalWeight", TotalEventWeight)

    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 



    # ErrorConditions
    # Undetectable (too low) event Nantennas = -5
    # Undetected (no trigger) candidate Nantennas = -4
    # AntennaP2PInfoNotFound Nantennas = -3 (means interpolation failed miserably)
    # Nantennas = 0 (means interpolation didnt actually get any antennas interpolated, probabbly becouse they where all discarded for being outside of the starshape, or the 4 antennas that where going to be used on the interpol where blow threshold)
    # Plane Reco Failed PlaneNAntennas = -1 (again, probbably becouse there where not enough antennas to start)
    # Sphere Reco Failed SphereNAntennas = -1
    # ADF Reco Failed ADFNAntennas = -1


    # STATISTICAL NOTE:
    # Since the method for generating the events tries core positions untill it cant find any more, the resulting probabbility of detection for a given zenith, energy bin (if we consider all the events to have the same chance)
    # is a geometric distribution, where we try k times until we succeed.
    # the formula for estimating the detection probability (which is also maximum likelyhood)is p = n/Sum(k) where n is the size . Variance of the distribution is 1-p/p^2
    # all this works if k is considered unlimited, that is, we will alway try until we get a success, and the probabbility of success is different from 0. 
    # in cases were we truncate the series prematurelly, or when we know from fact that the chance is 0 (The too low events) we are violating this.
    # The case where in a bin there are events that can be detected, and events that cannot, is an extreme case where is obvious that the probability is not equal for all events.
    # 
    # In this study, we have several steps that affect the detection chance
    # 0) Probability of an event to be undetectable (we postulate the existane of two populations. To ensure the number of events we are sampling is enough, check the distribution of the number of tries)
    # 1) Probabbility of Being a candidate (estimated with the geometric distro)
    # 1) The probabbility of a candidate to form a T2
    # 2) Probabbility of a T2 to be reconstructed (estimated as the succsessfull recos/total number of candidates) (where successfull can mean also successfull and good )
    # The final probability of detection is the product of the three? only if these probabilities are independent

    # Making the pandas dataframe
    df = pd.DataFrame(np.stack(( EventWeight , NAntennas , NTriggered , Energy , Zenith , Azimuth , PlaneNAntennas , PlaneZenith , PlaneAzimuth , PlaneChiSignif , PlaneChi2 , SphereNAntennas , SphereChiSignif , SphereChi2 , ADFNAntennas , ADFZenith , ADFAzimuth , ADFWidth , ADFAmplitude , ADFChiSignif , ADFChi2 ),axis = 1), 
        columns = ["EventWeight","NAntennas","NTriggered","Energy","Zenith","Azimuth","PlaneNAntennas","PlaneZenith","PlaneAzimuth","PlaneChiSignif","PlaneChi2","SphereNAntennas","SphereChiSignif","SphereChi2","ADFNAntennas","ADFZenith","ADFAzimuth","ADFWidth","ADFAmplitude","ADFChiSignif","ADFChi2"])

    
    # This is to add the series separating the points (im not using it for now, but it is recomended for performance)
    # pointsdf = pd.DataFrame(np.concatenate(( XmaxPosition , CorePosition , SphereSource),axis = 1),
    #       columns = ["XXmax","YXmax","ZXmax","XCore","YCore","ZCore","XSource","YSource","ZSource"])

    # This is to add the series with thepositions one by one, and in this cas keep the points in a list [x,y,z]
    df["XmaxPosition"] = XmaxPosition.tolist()
    df["CorePosition"] = CorePosition.tolist()
    # Accessing the first element of an array, for one condition
    # df[df["Energy"]>3.9]["CorePosition"].str[0]

    # Lets compute the pointing errors
    def ComputeAngularDistance(azim_r, zen_r, azim_s, zen_s) :
        '''
        Toi même tu sais
        '''
        # TODO: handle errors
        azim_diff = azim_r - azim_s
        return 180./np.pi * np.arccos(np.cos(zen_r*np.pi/180)*np.cos(zen_s*np.pi/180) + np.cos(azim_diff*np.pi/180) * np.sin(zen_s*np.pi/180) * np.sin(zen_r*np.pi/180))

    df["PlaneError"] = ComputeAngularDistance(PlaneAzimuth, PlaneZenith, Azimuth, Zenith)
    df["ADFError"] = ComputeAngularDistance(ADFAzimuth, ADFZenith, Azimuth, Zenith)
    # and the distance 2 core.
    df["Core2Center"] = np.sqrt(np.power(CorePosition[:,0]-Xc,2)+np.power(CorePosition[:,1]-Yc,2))

    df.to_csv('df_all_test.csv')

if False:
    read_files_and_save_csv()

df_bin = pd.read_csv('energ_zenith_bin.csv')
EnergyBins = df_bin["Energy"].unique()
EnergyBins.sort()
ZenithBins = df_bin["Zenith"].unique()
ZenithBins.sort()
ZenithBins = ZenithBins[::-1]  

df = pd.read_csv('df_all_test.csv')
loopcount = df.shape[0]
# Now, i have to incorporate the information of events that despite having the posibility in theory of producing triggers, they didnt.
# To consider an event to be a candidate, we ask that: 
# at least 4 antennas in the starshape are inside the array area AND above threshold (in E and V) 
# that the farthest away triggered antenna is on ring 2 (so, antenna 8 or more)) 
# (this avoids the problem of having the only triggered antennas to be on the "inner" circle, which on very close showers/very broad starshapes gives a very crappy intrpolation, makin all the inside of the inner circle basically constant.)
# If a starshape event failed to produce a viable event candidate after a number of tries (1000), it produces a file ending in .notrigger.txt, with all the tested core positions inside.
# but if the starshape had all the antennas below threshold,it will bever be detected, and a "toolow.txt" file is produced 

# Candidate events (considered as "detected" events below) are passed through the interpolation, but this is no guarantee that enough interpolation signals will be above threshold to form a T2 and then be granted a reconstruction 
# or that the interpolation will be successfull


# To see how well this candidate condition performs, lets check the Average Weight per Energy/Zenith bin (in a very panda power way!)
# limit plot to eight columns
MaxNtries = np.max(df["EventWeight"])-1
columns = 5
Bins = ZenithBins
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
    rows = rows + 1
# plot to see how is the average weight distributed, to detect problems
#df[(df["EventWeight"]>0)].groupby(["Energy","Zenith"])["EventWeight"].mean().unstack().plot(logx = True,subplots = True,sharex = True,sharey = True,layout = (rows,columns)) 
#plt.show()

# plot to see the histograms, and detect problems
#df[df["EventWeight"]>0]["EventWeight"].hist(by = df["Energy"],bins = int(1000/5),range = (0,1000),sharex = True)
#plt.show()

grouped = df.groupby(["Energy","Zenith"])
TotalEvents = []
EnergyCheck = []
ZenithCheck = []
Undetectable = []
Undetected = []
UndetectionChance = []
Detected = []
DetectedWeight = []
InterpolFailed = []
DetectionChance = []
T2Events = []
T2Chance = []
T2Efficiency = []
T2EffectiveArea = []
RecoFailed = []
ReconstructionChance = []
Labels = []
for label,gd in grouped:
    # 
    # gd = gd[gd["Core2Center"]<siteradius]
    # 
    print("# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # ")
    print("Analyzing Energy",label[0],"Zenith",180-label[1])
    print("# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # ")
    # Total Number of events
    Labels.append(label)
    TotalEvents.append(gd["NAntennas"].count())
    Events = TotalEvents[-1]
    EnergyCheck.append(label[0])
    ZenithCheck.append(label[1])
    if(Events<1):
        print("No Events On this bin")
        Undetectable.append(0)
        Undetected.append(0)
        UndetectionChance.append(1)
        Detected.append(0)
        DetectedWeight.append(0)
        InterpolFailed.append(0)
        DetectionChance.append(0)
        T2Events.append(0)
        T2Chance.append(0)
        T2Efficiency.append(0)
        T2EffectiveArea.append(0)
        RecoFailed.append(0)
        ReconstructionChance.append(0) 
        continue
 # 
    print("Bin Events:",Events,"(",100*Events/loopcount,"% of all events)")
    Undetectable.append(gd[gd["NAntennas"]  == -5]["NAntennas"].count())
    print("of wich Undetectable (toolow):",Undetectable[-1],"(",100*Undetectable[-1]/Events,"%)")
    Undetected.append(gd[gd["NAntennas"]  == -4]["NAntennas"].count())
    UndetectionChance.append((Undetectable[-1]+Undetected[-1])/Events)
    print("and Undetected (notrigger):",Undetected[-1],"(",100*Undetected[-1]/Events,"%), establishing an undetection chance of",100*UndetectionChance[-1],"%")
    Detected.append(Events-Undetectable[-1]-Undetected[-1])
    Detected2 = gd[(gd["NAntennas"]!= -4)&(gd["NAntennas"]!= -5)]["EventWeight"].count() 
    InterpolFailed.append(gd[(gd["NAntennas"]!= -4)&(gd["NAntennas"]!= -5)&(gd["NAntennas"]  == -3)]["EventWeight"].count())
    if(Detected[-1]!= Detected2):
        print("There is a mismatch between detected and detected2, check why",InterpolFailed[-1],Detected[-1],Detected2)
        DetectedWeight.append(0)
        InterpolFailed.append(0)
        DetectionChance.append(0)
        T2Events.append(0)
        T2Chance.append(0)
        T2Efficiency.append(0)
        T2EffectiveArea.append(0)
        RecoFailed.append(0) 
        ReconstructionChance.append(0) 
        continue
 #  
    if(InterpolFailed[-1]>0): 
        print("and the Interpolation Failed on",InterpolFailed[-1],"which we ignore")
        Events = Events-InterpolFailed[-1]
        Detected[-1] = Events-Undetectable[-1]-Undetected[-1]
 # 
    DetectedWeight.append(gd[(gd["NAntennas"]!= -4)&(gd["NAntennas"]!= -5)&(gd["NAntennas"]!= -3)]["EventWeight"].sum())
    if(DetectedWeight[-1]<1):
        print("No Events Detected on this bin")
        DetectionChance.append(0)
        T2Events.append(0)
        T2Chance.append(0)
        T2Efficiency.append(0)
        T2EffectiveArea.append(0)
        RecoFailed.append(0)
        ReconstructionChance.append(0)
        continue 
        # 
    DetectionChance.append(Detected[-1]/DetectedWeight[-1])
    print("This leaves",Detected[-1],"Detected events, that represent",DetectedWeight[-1],"trials, impling tha Detectable events have a ",100*DetectionChance[-1],"% chance of being actually detected, depending on where they landed")
    print("these ",Detected[-1],"detected events, represent",Detected[-1]/Events,"% of the total number of simulated events")
    # 
    T2Events.append(gd[(gd["NAntennas"]!= -4)&(gd["NAntennas"]!= -5)&(gd["NAntennas"]!= -3)&(gd["NTriggered"]>= TriggerN)]["EventWeight"].count())
    T2Chance.append(T2Events[-1]/Detected[-1])
    NoT2Events = Detected[-1]-T2Events[-1]
    print("From these,",T2Events[-1],"formed a T2 trigger, giving a ",100*T2Chance[-1],"% Chance of T2 detection")
    T2Efficiency.append(T2Chance[-1]*DetectionChance[-1])
    step = 5750/np.cos((180-label[1])*np.pi/180.0) # step used on event generation (random points in an exagon of side "step")
    SampledArea = 3.0*np.sqrt(3)*step*step/2.0 # area of the exagon of side step
    T2EffectiveArea.append(SampledArea*T2Efficiency[-1])
    print("This gives a T2Efficiency of",100*T2Efficiency[-1],"% and T2EffectiveArea of",T2EffectiveArea[-1]/1E6,"[km2]. DetectorArea:",DetectorArea/1E6)
    if(T2Events[-1]  == 0):
        RecoFailed.append(0)
        ReconstructionChance.append(0)
   
    PlaneRecoFailed = gd[ (gd["PlaneNAntennas"]  == -1) & (gd["NTriggered"]>= TriggerN) ]["NAntennas"].count()
    # print("Of the Detected,",PlaneRecoFailed,"Failed the PlaneReco",100*PlaneRecoFailed/Events,"% of all",100*PlaneRecoFailed/Detected[-1],"%of the Detected")
    SphereRecoFailed = gd[ (gd["SphereNAntennas"]  == -1) & (gd["NTriggered"]>= TriggerN) ]["NAntennas"].count()
    # print("Of the Detected,",SphereRecoFailed,"Failed the SphereReco",100*SphereRecoFailed/Events,"% of all",100*SphereRecoFailed/Detected[-1],"% of the Detected")
    ADFRecoFailed = gd[(gd["ADFNAntennas"]  == -1)&(gd["NTriggered"]>= TriggerN)]["NAntennas"].count()
    # print("Of the Detec,",ADFRecoFailed,"Failed the ADFReco",100*ADFRecoFailed/Events,"% of all",100*ADFRecoFailed/Detected[-1],"% of the Detected")
    RecoFailed.append(gd[((gd["PlaneNAntennas"]  == -1)|(gd["SphereNAntennas"]  == -1)|(gd["ADFNAntennas"]  == -1))&(gd["NTriggered"]>= TriggerN)]["NAntennas"].count())
    print("Of the Detected,",RecoFailed[-1],"Failed the Reco",100*RecoFailed[-1]/T2Events[-1],"% of all",100*RecoFailed[-1]/Detected[-1],"% of the Detected")
    Reconstructed = T2Events[-1]-RecoFailed[-1]
    ReconstructionChance.append(Reconstructed/T2Events[-1])
    print("so we are left with",Reconstructed,"Reconstructed events, implying that Detected have a ",100*Reconstructed/T2Events[-1],"of being reconstructed successfully" ) 
    ReconstructionEfficiency = ReconstructionChance[-1]*T2Chance[-1]*DetectionChance[-1]
    # sumcheck
    SumCheck = Undetectable[-1]+Undetected[-1]+InterpolFailed[-1]+NoT2Events+RecoFailed[-1]+Reconstructed 
    if(SumCheck!= TotalEvents[-1]):
        print("WARNING SOMETHING DOES NOT ADD UP, CHECK THIS")
    

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
# making the Bining for the pcolormesh plots
Labels = np.array(Labels) # wq have a 0,0 group that i will ignore
EnergyCheck = np.array(EnergyCheck)
ZenithCheck = np.array(ZenithCheck)

Energies = np.unique(Labels[:,0])
LogEnergies = np.log10(Energies*1E18)
DeltaLogEnergies = np.abs(np.average((LogEnergies[:-1]-LogEnergies[1:])[1:]))/2
LowEdgeLogEnergies = LogEnergies[1:]-DeltaLogEnergies
LowEdgeLogEnergies = np.append(LowEdgeLogEnergies,LogEnergies[-1]+DeltaLogEnergies)

Zeniths = 180.0-np.unique(Labels[:,1])
LogSecZeniths = np.log10(1.0/np.cos(np.deg2rad(Zeniths)))
DeltaLogSecZeniths = np.abs(np.average((LogSecZeniths[:-1]-LogSecZeniths[1:])[1:]))/2
LowEdgeLogSecZeniths = LogSecZeniths[1:]+DeltaLogSecZeniths
LowEdgeLogSecZeniths = np.append(LowEdgeLogSecZeniths,LogSecZeniths[-1]-DeltaLogSecZeniths)

LowEdgeZeniths = np.rad2deg(np.arccos(np.power(10,-LowEdgeLogSecZeniths)))
SolidAngles = 2*np.pi*np.cos(np.deg2rad(LowEdgeZeniths[1:]))-np.cos(np.deg2rad(LowEdgeZeniths[:-1]))

X, Y = np.meshgrid(LowEdgeLogEnergies, LowEdgeLogSecZeniths)

MeshZeniths = (180-Labels[:,1])
MeshCos = np.cos(np.deg2rad(MeshZeniths))
MeshLogSecZeniths = np.log10(1.0/np.cos(np.deg2rad(MeshZeniths)))
MeshSolidAngles = 2*np.pi*(np.power(10,-(MeshLogSecZeniths-DeltaLogSecZeniths))-np.power(10,-(MeshLogSecZeniths+DeltaLogSecZeniths)))

width = 10
height = width/1.618
plt.rc('font', family = 'serif', size = 15)
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 

def PhaseSpacePlot(ax,PlotTitle,X,Y,Z):
 tmp = ax.pcolormesh(X, Y, Z.T,edgecolors = 'k', linewidth = 2)
 tmp = fig100.colorbar(tmp, ax = ax,pad = 0.15)
 # tmp = ax.set_aspect('equal')
 tmp = ax.set_title(PlotTitle)
 tmp = ax.set_xlabel("$Log_{10}(Energy[eV])$")
 tmp = ax.set_ylabel("$Log_{10}(sec(Zenith))$")
 # for secondary axis
 lambdaZenith = lambda LSZ: np.rad2deg(np.arccos(np.power(10,-LSZ)))
 lambdaLSZ = lambda Zenith: np.log10(1.0/np.cos(np.deg2rad(Zenith)))
 ax2 = ax.secondary_yaxis('right',functions = (lambdaZenith,lambdaLSZ))
 ax2.set_yticks(Zeniths[1:])
 ax2.set_yticklabels([str(round(float(label), 1)) for label in Zeniths[1:]]) 
 tmp = ax2.set_ylabel("$Zenith [deg]$")
 # 
 
# TotalEvents
# Z = np.reshape(TotalEvents[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(TotalEvents[0:], (len(EnergyBins)-1, len(ZenithBins)-1))
fig100 = plt.figure(100, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig100.set_tight_layout(True)
ax = fig100.add_subplot(111)

PhaseSpacePlot(ax,"Total Events",X,Y,Z)

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # -
# EnergyCheck
fig101 = plt.figure(101, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
Z = np.reshape(np.log10(EnergyCheck[0:]),(len(EnergyBins)-1,len(ZenithBins)-1))
# Z = np.reshape(np.log10(EnergyCheck[1:-1]),(len(EnergyBins)-1,len(ZenithBins)-2))
fig101 = plt.figure(101, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig101.set_tight_layout(True)
ax = fig101.add_subplot(221)

PhaseSpacePlot(ax,"$Log_{10}(Energy[eV])$ Check",X,Y,Z)

# ZenithCheck
# Z = np.reshape(180-ZenithCheck[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(180-ZenithCheck[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig101.add_subplot(222)
PhaseSpacePlot(ax,"Zenith Check",X,Y,Z)

# InterpolFailed
# Z = np.reshape(InterpolFailed[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(InterpolFailed[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig101.add_subplot(223)
PhaseSpacePlot(ax,"InterpolFailed",X,Y,Z)

# Z = np.reshape(RecoFailed[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(RecoFailed[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig101.add_subplot(224)
PhaseSpacePlot(ax,"RecoFailed",X,Y,Z)

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Undetectable
fig102 = plt.figure(102, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig102.set_tight_layout(True)
ax = fig102.add_subplot(221)

# Z = np.reshape(Undetectable[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(Undetectable[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
PhaseSpacePlot(ax,"Number Undetectable",X,Y,Z)

# Undetected
# Z = np.reshape(Undetected[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(Undetected[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig102.add_subplot(222)
PhaseSpacePlot(ax,"Number Undetected",X,Y,Z)

Z1 = np.reshape(Undetectable[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z2 = np.reshape(Undetected[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = Z1+Z2
ax = fig102.add_subplot(223)
PhaseSpacePlot(ax,"Undetected + Undetectable",X,Y,Z)

# UndetectionChance
# Z = np.reshape(UndetectionChance[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = 100*np.reshape(UndetectionChance[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig102.add_subplot(224)
PhaseSpacePlot(ax,"Undetection Chance [%]",X,Y,Z)


# T2Chance
# Z = np.reshape(T2Chance[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
# Z = np.reshape(T2Chance[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
# Z = 1-Z
# ax = fig102.add_subplot(224)
# PhaseSpacePlot(ax,"Untrigger (Candidate -> T2) %",X,Y,Z)


# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Detected
fig103 = plt.figure(103, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig103.set_tight_layout(True)

# Candidate
# Z = np.reshape(Detected[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(Detected[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig103.add_subplot(131)
PhaseSpacePlot(ax,"Candidated",X,Y,Z)

# DetectedWeight
# Z = np.reshape(DetectedWeight[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(DetectedWeight[0:],(len(EnergyBins)-1,len(ZenithBins)-1))

ax = fig103.add_subplot(132)
PhaseSpacePlot(ax,"Candidated Weight",X,Y,Z)

# DetectionChance = []
# Z = np.reshape(DetectionChance[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(DetectionChance[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig103.add_subplot(133)
PhaseSpacePlot(ax,"Candidated Chance",X,Y,Z)

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 3
# T2Events
fig104 = plt.figure(104, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig104.set_tight_layout(True)

# DetectionChance = []
# Z = np.reshape(DetectionChance[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = 100*np.reshape(DetectionChance[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig104.add_subplot(221)
PhaseSpacePlot(ax,"Candidated Chance",X,Y,Z)

# T2Chance
# Z = np.reshape(T2Chance[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = 100*np.reshape(T2Chance[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig104.add_subplot(222)
PhaseSpacePlot(ax,"T2Chance (Candidate -> T2) %",X,Y,Z)

# Z = np.reshape(T2Efficiency[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = 100*np.reshape(T2Efficiency[0:],(len(EnergyBins)-1,len(ZenithBins)-1))

ax = fig104.add_subplot(223)
PhaseSpacePlot(ax,"T2Acceptance: (Event -> T2) %",X,Y,Z)


# Z = np.reshape(T2Events[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(T2Events[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig104.add_subplot(224)
PhaseSpacePlot(ax,"Number of T2Events",X,Y,Z)


# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 3

fig105 = plt.figure(105, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig105.set_tight_layout(True)


# Cosine
# Z = np.reshape(MeshCos[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(MeshCos[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig105.add_subplot(221)
PhaseSpacePlot(ax,"Projection Cosine",X,Y,Z)

# SolidAngle
# Z = np.reshape(MeshSolidAngles[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(MeshSolidAngles[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig105.add_subplot(222)
PhaseSpacePlot(ax,"Solid Angles",X,Y,Z)


T2EffectiveAreaCos = T2EffectiveArea*MeshCos
# Z = np.reshape(T2EffectiveAreaCos[1:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
Z = np.reshape(T2EffectiveAreaCos[0:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
ax = fig105.add_subplot(223)
PhaseSpacePlot(ax,"T2Area x Direction Projection $[km^2]$ ",X,Y,Z)


T2EffectiveAreaCosSr = T2EffectiveArea*MeshCos*MeshSolidAngles
# Z = np.reshape(T2EffectiveAreaCos[1:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
Z = np.reshape(T2EffectiveAreaCosSr[0:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
ax = fig105.add_subplot(224)
PhaseSpacePlot(ax,"Aperture $[km^2.sr]$",X,Y,Z)

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 3
# Money Plots
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 3# # 

fig106 = plt.figure(106, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig106.set_tight_layout(True)

# Z = np.reshape(T2Efficiency[1:],(len(EnergyBins)-1,len(ZenithBins)-1))
Z = np.reshape(T2Efficiency[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
ax = fig106.add_subplot(221)
PhaseSpacePlot(ax,"T2Acceptance: (Event -> T2) %",X,Y,Z)

# T2EffectiveArea
# Z = np.reshape(T2EffectiveArea[1:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
Z = np.reshape(T2EffectiveArea[0:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
ax = fig106.add_subplot(222)
PhaseSpacePlot(ax,"T2EffectiveArea $[km^2]$",X,Y,Z)

T2EffectiveAreaCos = T2EffectiveArea*MeshCos
Z = np.reshape(T2EffectiveAreaCos[0:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
ax = fig106.add_subplot(223)
PhaseSpacePlot(ax,"T2Area x Direction Projection $[km^2]$ ",X,Y,Z)

T2EffectiveAreaCosSr = T2EffectiveArea*MeshCos*MeshSolidAngles
Z = np.reshape(T2EffectiveAreaCosSr[0:],(len(EnergyBins)-1,len(ZenithBins)-1))/1E6
ax = fig106.add_subplot(224)
PhaseSpacePlot(ax,"Aperture $[km^2.sr]$",X,Y,Z)


Z = np.reshape(T2EffectiveAreaCosSr[0:],(len(Energies),len(Zeniths)))/1E6
ApertureZenith = Z.sum(axis = 0)
ApertureEnergy = Z.sum(axis = 1)

fig107 = plt.figure(107, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig107.set_tight_layout(True)
ax = fig107.add_subplot(121)
tmp = ax.plot(Zeniths[0:],ApertureZenith)
tmp = ax.set_xlabel("Zenith $[deg]$")
tmp = ax.set_ylabel("Aperture $[km^2.sr]$")
tmp = ax.set_ylim(0,140)
ax = fig107.add_subplot(122)
tmp = ax.plot(Energies[0:],ApertureEnergy)
tmp = ax.set_xlabel("Energy $[EeV]$")
tmp = ax.set_xscale('log')
tmp = ax.set_ylabel("Aperture $[km^2.sr]$")
tmp = ax.set_ylim(0,140)

# Mirar http://arpg-serv.ing2.uniroma1.it/patera/didattica/fis_mod/trasp_riv/Glossario/node3.html para ver como tratar los errores

 
def TALEDiferentialFlux(Energy):# eV⁻¹.m⁻²sr⁻¹s⁻¹
 # expects input in GeV
 # spectrum data from TALE 2018 https://arxiv.org/pdf/1803.01288.pdf (note that this was updated on ICRC2021 here https://pos.sissa.it/395/347/pdf work to be done!)
 # break point at log10E 16.22 slope -2.92 up to 17.04
 # at 16.2 the diferential flux is E³J(E) = 2.978 ±0.030 + 0.673 - 0.526 x10E24 eV²m⁻²sr⁻¹s⁻¹]
 # at 17 is 3.409 ± 0.113 + 0.832 - 0.482
 # J = A*E^gamma
 # double flux = A*pow(1E17,gamma);
 # double fluxE3 = flux*pow(1E17,3);
 # printf("A %g flux17 %g flux17*E³ %g\n",A,flux,fluxE3);
 # flux = A*pow(pow(10,15.5),gamma);
 # fluxE3 = flux*pow(pow(10,15.5),3);
 # printf("A %g flux15.5 %g flux15.5*E³ %g\n",A,flux,fluxE3);
 gamma2 = -2.92
 if((Energy>= np.power(10,8.05)) and (Energy<np.power(10,9.7))):
  gamma3 = -3.19
  J3E3 = 3.504*1E24
  E3 = np.power(10,17.05)
  A3 = (J3E3/np.power(E3,3))/np.power(E3,gamma3) 
  Flux3 = A3*np.power(Energy*1E9,-3.19)
  return Flux3
 elif((Energy>= pow(10,7.2)) and (Energy<pow(10,8.05))):  
  J2E3 = 3.504*1E24
  E2 = np.power(10,17.05)
  A2 = (J2E3/np.power(E2,3))/np.power(E2,gamma2)
  Flux2 = A2*np.power(Energy*1E9,gamma2) 
  return Flux2
 # break point at log10E 16.22 slope -2.92 up to 17.04
 # at 16.2 the diferential flux is E³J(E) = 2.978 ±0.030 + 0.673 - 0.526 x10E24 eV²m⁻²sr⁻¹s⁻¹]
 # at 15.7 3.405 ± 0.022 + 0.413 - 0.413
 elif( (Energy>= np.power(10,6.7)) and (Energy<np.power(10,7.2))):
  gamma1 = -3.12
  J1E3 = 3.405*1E24
  E1 = np.power(10,15.7)
  A1 = (J1E3/np.power(E1,3))/np.power(E1,gamma1) 
  Flux1 = A1*np.power(Energy*1E9,gamma1)
  return Flux1
 # //and then i assume a flux with gamma0 = gamma2;
 elif((Energy>np.power(10,6)) and (Energy < np.power(10,6.7))):
  gamma0 = gamma2
  J0E3 = 3.405*1E24
  E0 = np.power(10,15.7)
  A0 = (J0E3/np.power(E0,3))/np.power(E0,gamma0)
  Flux0 = A0*np.power(Energy*1E9,gamma0)
  return Flux0
 else:
  print("flux is not defined in this energy range",np.log10(Energy))
  return 0
 
xvalues = np.power(10,np.arange(6.1,9.4,0.1))
yvalues = np.zeros(len(xvalues))
for i,Energy in enumerate(xvalues):
 yvalues[i] = TALEDiferentialFlux(Energy)*np.power(Energy*1E9,3)

plt.scatter(xvalues,yvalues)
plt.loglog()
plt.show() 


Flux = np.zeros(len(EnergyCheck))
EventRate = np.zeros(len(EnergyCheck))
EventsDay = np.zeros(len(EnergyCheck))

for i,Energy in enumerate(EnergyCheck):
 Flux[i] = TALEDiferentialFlux(Energy*1E9) # eV⁻¹.m⁻²sr⁻¹s⁻¹
 DeltaE = (np.power(10,np.log10(Energy)+DeltaLogEnergies) -np.power(10,np.log10(Energy)-DeltaLogEnergies))*1E18 # eV
 EventRate[i] = Flux[i]*DeltaE*T2EffectiveAreaCosSr[i] # s⁻¹
 EventsDay[i] = EventRate[i]*24*60*60 # day⁻¹


Z = np.reshape(EventsDay[0:],(len(EnergyBins)-1,len(ZenithBins)-1))
fig108 = plt.figure(108, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig108.set_tight_layout(True)
ax = fig108.add_subplot(111)
PhaseSpacePlot(ax,"Event Rate: (Events/day)",X,Y,Z) 
 

Z = np.reshape(EventsDay[0:],(len(Energies),len(Zeniths)))
EventRateZenith = Z.sum(axis = 0)
EventRateEnergy = Z.sum(axis = 1)
print("Check TotalEventRate coincide ",EventRateZenith.sum(),EventRateEnergy.sum())

fig107 = plt.figure(107, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
fig107.set_tight_layout(True)
ax = fig107.add_subplot(121)
tmp = ax.plot(Zeniths[0:],EventRateZenith)
tmp = ax.set_xlabel("Zenith $[deg]$")
tmp = ax.set_ylabel("Event Rate $[day^{-1}]$")
tmp = ax.set_yscale('log')
# tmp = ax.set_ylim(0,140)
ax = fig107.add_subplot(122)
tmp = ax.plot(Energies[0:],EventRateEnergy)
tmp = ax.set_xlabel("Energy $[EeV]$")
tmp = ax.set_xscale('log')
tmp = ax.set_yscale('log')
tmp = ax.set_ylabel("Event Rate $[day^{-1}]$")
# tmp = ax.set_ylim(0,140)



 
 
# OtherPlots

df[((df["PlaneNAntennas"]  == -1)|(df["SphereNAntennas"]  == -1)|(df["ADFNAntennas"]  == -1))&(df["NTriggered"]>= TriggerN)]


df[(df["NTriggered"]>5)&(df["Zenith"]<113)]["Core2Center"].hist(by = 180-df["Zenith"],bins = 100,sharex = True)

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
# 2D Plot of core positions (for checking)
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 

fig101 = plt.figure(101, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
ax1 = fig101.add_subplot(121)
xvalues = CorePosition[:,0]
yvalues = CorePosition[:,1]
zvalues = Azimuth
im = ax1.scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(zvalues),vmax = np.max(zvalues), s = 4,alpha = 0.5)
bar = fig101.colorbar(im,ax = ax1)
bar.set_label('Azimuth', rotation = 270, labelpad = 17)
tmp = ax1.set_ylabel('Westing[m]')
tmp = ax1.set_xlabel("Northing[m]")
circle2 = plt.Circle((Xc, Yc), siteradius, color = 'green',alpha = 0.25)
ax1.add_patch(circle2)

ax2 = fig101.add_subplot(122)
xvalues = CorePosition[:,0]
yvalues = CorePosition[:,1]
zvalues = 180-Zenith
im = ax2.scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(zvalues),vmax = np.max(zvalues), s = 4,alpha = 0.5)
bar = fig101.colorbar(im,ax = ax2)
bar.set_label('Zenith', rotation = 270, labelpad = 17)
tmp = ax2.set_ylabel('Westing[m]')
tmp = ax2.set_xlabel("Northing[m]")
circle2 = plt.Circle((Xc, Yc), siteradius, color = 'green',alpha = 0.25)
ax2.add_patch(circle2)

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
# 3D Plot
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
rearth = 6370949.0

from mpl_toolkits.mplot3d import Axes3D
fig102 = plt.figure(102, figsize = (width,height), facecolor = 'w', edgecolor = 'k')
ax1 = fig102.add_subplot(121,projection = '3d')
xvalues = XmaxPosition[:,0]
yvalues = XmaxPosition[:,1]
zvalues = XmaxPosition[:,2]

R02 = np.power(xvalues,2)+np.power(yvalues,2)
altitude = (np.sqrt(np.power(zvalues+rearth,2) + R02 ) - rearth) # altitude of emission, in km

colors = Azimuth
im = ax1.scatter(xvalues,yvalues,altitude,c = colors,cmap = plt.cm.jet,vmin = np.min(colors),vmax = np.max(colors), alpha = 0.5)
bar = fig102.colorbar(im,ax = ax1)
bar.set_label('Azimuth', rotation = 270, labelpad = 17)
tmp = ax1.set_ylabel('Westing[m]')
tmp = ax1.set_xlabel("Northing[m]")

ax2 = fig102.add_subplot(122,projection = '3d')
xvalues = XmaxPosition[:,0]
yvalues = XmaxPosition[:,1]
zvalues = XmaxPosition[:,2]
colors = 180-Zenith
im = ax2.scatter(xvalues,yvalues,altitude,c = colors,cmap = plt.cm.jet,vmin = np.min(colors),vmax = np.max(colors),alpha = 0.5)
bar = fig102.colorbar(im,ax = ax2)
bar.set_label('Zentih', rotation = 270, labelpad = 17)
tmp = ax2.set_ylabel('Westing[m]')
tmp = ax2.set_xlabel("Northing[m]")





# this will plot two variables 
ax = df.plot("ADFError","PlaneError",kind = "scatter")
# this will plot two variables with a condition
ax1 = df[df["Core2Center"]<siteradius].plot("ADFError","PlaneError",kind = "scatter")
ax2 = df[df["Core2Center"]>siteradius].plot("ADFError","PlaneError",kind = "scatter")
plt.show()
# Two conditions
# ax1 = df[(df["Core2Center"]<siteradius) & (df["Energy"]>3.9)].plot("ADFError","PlaneError",kind = "scatter")



# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 

# So now, lets try to plot the core positions of the successfully reconstructed events
ErrorThreshold = 1

# Getting unique memebers of a Series
EnergyBins = df["Energy"].unique()
EnergyBins.sort()
ZenithBins = df["Zenith"].unique()
ZenithBins.sort()
ZenithBins = ZenithBins[::-1] # reverse the order

# Getting unique combination of a combination of series
# uniquebins = df.drop_duplicates(["Zenith","Energy"])
# ZenithAzimuthBins = uniquebins[["Energy","Zenith"]


# limit plot to eight columns
columns = 4
Bins = ZenithBins[1:-1]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
ErrorThreshold = 1000 
for n, Bin in enumerate(Bins): 
 xvalues = df[(df["PlaneError"]<ErrorThreshold)&(df["Zenith"]  == Bin)&(df["NTriggered"]>4)&(df["PlaneError"]<1)&(df["ADFError"]<2.5)]["CorePosition"].str[0]
 yvalues = df[(df["PlaneError"]<ErrorThreshold)&(df["Zenith"]  == Bin)&(df["NTriggered"]>4)&(df["PlaneError"]<1)&(df["ADFError"]<2.5)]["CorePosition"].str[1]
 zvalues = df[(df["PlaneError"]<ErrorThreshold)&(df["Zenith"]  == Bin)&(df["NTriggered"]>4)&(df["PlaneError"]<1)&(df["ADFError"]<2.5)]["ADFError"]
 # xvalues = df[df["Zenith"]  == Bin]["CorePosition"].str[0]
 # yvalues = df[df["Zenith"]  == Bin]["CorePosition"].str[1]
 # zvalues = df[df["Zenith"]  == Bin]["Azimuth"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(0),vmax = np.max(1), s = 4,alpha = 0.5)
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('Westing[m]')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("Northing[m]")
 tmp = ax1.flat[n].set_title("Zenith "+"%.1f" % (180 - Bin))
 circle2 = plt.Circle((Xc, Yc), siteradius, color = 'yellow',alpha = 0.2)
 tmp = ax1.flat[n].add_patch(circle2) 



# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #  

# Now, lets plot the angular error correlation for each energy bin

# limit plot to eight columns
columns = 4
Bins = EnergyBins[1:]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
fig1.set_tight_layout(True)  
for n, Bin in enumerate(Bins): 
 yvalues = df[(df["Energy"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneError"]
 xvalues = df[(df["Energy"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["ADFError"]
 zvalues = df[(df["Energy"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneNAntennas"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(TriggerN),vmax = np.max(TriggerN*6), s = 4,alpha = 0.5)
 # bar = fig1.colorbar(im,ax = ax1.flat[n])
 # bar.set_label('Nantennas', rotation = 270, labelpad = 17)
 tmp = ax1.flat[n].set_ylim(0,1)
 tmp = ax1.flat[n].set_xlim(0,2.5) 
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('PlaneError[deg]')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("ADFError[deg]")
 tmp = ax1.flat[n].set_title("Energy "+"%.1f" % (18.0+np.log10(Bin)))


# and each zenith bin (color NAnt)

columns = 4
Bins = ZenithBins[1:-1]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
 
for n, Bin in enumerate(Bins): 
 yvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneError"]
 xvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["ADFError"]
 zvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneNAntennas"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(TriggerN),vmax = np.max(TriggerN*6), s = 4,alpha = 0.5)
 # bar = fig1.colorbar(im,ax = ax1.flat[n])
 # bar.set_label('Nantennas', rotation = 270, labelpad = 17)
 tmp = ax1.flat[n].set_ylim(0,1)
 tmp = ax1.flat[n].set_xlim(0,2.5) 
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('PlaneError[deg]')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("ADFError[deg]")
 tmp = ax1.flat[n].set_title("Zenith"+"%.1f" % float(180.0-Bin)) 


# and each zenith bin (color Plane Error)

columns = 4
Bins = ZenithBins[1:-1]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

MinAntennasvsZenith = [ 5 ,  7,  8, 10, 15, 22, 24, 27, 30, 33, 35, 37]
#      [67.75, 71.61, 74.76, 77.35, 79.49, 81.26, 82.72, 83.94, 84.95, 85.8 , 86.5 , 87.08]


fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
  
for n, Bin in enumerate(Bins): 
 # yvalues = 18.0+np.log10(df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["Energy"])
 yvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["ADFError"]
 xvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneNAntennas"]
 zvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneError"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(0),vmax = np.max(0.5), s = 4,alpha = 0.5)
 bar = fig1.colorbar(im,ax = ax1.flat[n])
 bar.set_label('PlaneError', rotation = 270, labelpad = 17)
 # tmp = ax1.flat[n].set_ylim(0,1)
 # tmp = ax1.flat[n].set_xlim(0,2.5) 
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('ADFError')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("NAntennas")
 tmp = ax1.flat[n].axvline(x = MinAntennasvsZenith[n]) 
 tmp = ax1.flat[n].set_title("Zenith"+"%.1f" % float(180.0-Bin)) 



# And now with variable antenna cut

columns = 4
Bins = ZenithBins[1:-1]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
 
for n, Bin in enumerate(Bins): 
 print(180-Bin,MinAntennasvsZenith[n])
 yvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>= MinAntennasvsZenith[n])]["PlaneError"]
 xvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>= MinAntennasvsZenith[n])]["ADFError"]
 zvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>= MinAntennasvsZenith[n])]["PlaneNAntennas"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(TriggerN),vmax = np.max(TriggerN*6), s = 4,alpha = 0.5)
 # bar = fig1.colorbar(im,ax = ax1.flat[n])
 # bar.set_label('Nantennas', rotation = 270, labelpad = 17)
 tmp = ax1.flat[n].set_ylim(0,1)
 tmp = ax1.flat[n].set_xlim(0,2.5) 
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('PlaneError[deg]')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("ADFError[deg]")
 tmp = ax1.flat[n].set_title("Zenith "+"%.1f" % float(180.0-Bin)) 











# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


columns = 4
Bins = EnergyBins[1:]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
fig1.set_tight_layout(True) 
for n, Bin in enumerate(Bins): 
 xvalues = df[(df["Energy"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["ADFError"]
 yvalues = np.log10(1/np.cos(np.deg2rad(180-df[(df["Energy"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["Zenith"])))
 zvalues = df[(df["Energy"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneNAntennas"]
 # xvalues = df[df["Zenith"]  == Bin]["CorePosition"].str[0]
 # yvalues = df[df["Zenith"]  == Bin]["CorePosition"].str[1]
 # zvalues = df[df["Zenith"]  == Bin]["Azimuth"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(TriggerN),vmax = np.max(TriggerN*3), s = 4,alpha = 0.5)
 bar = fig1.colorbar(im,ax = ax1.flat[n])
 bar.set_label('Nantennas', rotation = 270, labelpad = 17)
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('log10(sec(zenith))')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("PlaneError[deg]")
 tmp = ax1.flat[n].set_title("Energy "+"%.1f" % (18.0+np.log10(Bin)))
 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #  
 
columns = 4
Bins = ZenithBins[:-1]
quotient, remainder = divmod(len(Bins), columns)
rows = quotient
if remainder:
 rows = rows + 1

fig1, ax1 = plt.subplots(rows, columns, sharex = 'all', sharey = 'all',figsize = (width,height)) 
# doing ax1.flat[n] i can go through the plots
 
for n, Bin in enumerate(Bins): 
 xvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneError"]
 yvalues = np.log10(df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["Energy"])+18
 zvalues = df[(df["Zenith"]  == Bin)&(df["PlaneNAntennas"]>TriggerN)]["PlaneNAntennas"]
 # xvalues = df[df["Zenith"]  == Bin]["CorePosition"].str[0]
 # yvalues = df[df["Zenith"]  == Bin]["CorePosition"].str[1]
 # zvalues = df[df["Zenith"]  == Bin]["Azimuth"]
 im = ax1.flat[n].scatter(xvalues,yvalues,c = zvalues,cmap = plt.cm.jet,vmin = np.min(TriggerN),vmax = np.max(TriggerN*3), s = 4,alpha = 0.5)
 bar = fig1.colorbar(im,ax = ax1.flat[n])
 bar.set_label('Nantennas', rotation = 270, labelpad = 17)
 quotient, remainder = divmod(n, columns)
 if(remainder  == 0):
  tmp = ax1.flat[n].set_ylabel('Energy')
 if(n>= columns*(rows-1)): 
  tmp = ax1.flat[n].set_xlabel("PlaneError[deg]")
 tmp = ax1.flat[n].set_title("Zenith"+"%.1f" % float(180.0-Bin)) 
 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
# 
# NOW WE REPEAT BUT with core
# 
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  
 
 
