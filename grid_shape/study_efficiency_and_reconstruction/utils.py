
import pandas as pd
import os
import glob
import numpy as np
import hdf5fileinout as hdf5io


def read_files_and_save_csv(
    library_base_directory,
    output_file1,
    output_file2,
    TriggerThreshold=45,
    TriggerN=5,
    Xc=0,
    Yc=0
):

    EXT = "*.hdf5"
    all_hdf5_files = [
        file for path, subdir, files in os.walk(library_base_directory) for file in glob.glob(os.path.join(path, EXT))
    ][::1]

    EXT = "*.notrigger.txt"
    all_notrigger_files = [
        file for path, subdir, files in os.walk(library_base_directory) for file in glob.glob(os.path.join(path, EXT))
    ][::1]

    EXT = "*.toolow.txt"
    all_toolow_files = [
        file for path, subdir, files in os.walk(library_base_directory) for file in glob.glob(os.path.join(path, EXT))
    ][::1]

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

    df.to_csv(output_file1)

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
    df = pd.DataFrame(np.stack((EventWeight, NAntennas, NTriggered, Energy, Zenith, Azimuth, PlaneNAntennas, PlaneZenith, PlaneAzimuth, PlaneChiSignif, PlaneChi2, SphereNAntennas, SphereChiSignif, SphereChi2, ADFNAntennas, ADFZenith, ADFAzimuth, ADFWidth, ADFAmplitude, ADFChiSignif, ADFChi2), axis=1), 
        columns = ["EventWeight","NAntennas","NTriggered","Energy","Zenith","Azimuth","PlaneNAntennas","PlaneZenith","PlaneAzimuth","PlaneChiSignif","PlaneChi2","SphereNAntennas","SphereChiSignif","SphereChi2","ADFNAntennas","ADFZenith","ADFAzimuth","ADFWidth","ADFAmplitude","ADFChiSignif","ADFChi2"])

    
    # This is to add the series separating the points (im not using it for now, but it is recomended for performance)
    # pointsdf = pd.DataFrame(np.concatenate(( XmaxPosition , CorePosition , SphereSource),axis = 1),
    #       columns = ["XXmax","YXmax","ZXmax","XCore","YCore","ZCore","XSource","YSource","ZSource"])

    # This is to add the series with thepositions one by one, and in this cas keep the points in a list [x,y,z]
    df["XmaxPosition_x"] = XmaxPosition[:,0]
    df["XmaxPosition_y"] = XmaxPosition[:,1]
    df["XmaxPosition_z"] = XmaxPosition[:,2]

    df["CorePosition_x"] = CorePosition[:,0]
    df["CorePosition_y"] = CorePosition[:,1]
    df["CorePosition_z"] = CorePosition[:,2]
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

    df.to_csv(output_file2)

