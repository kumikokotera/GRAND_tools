'''Script to perform an interpolation between to electric field traces at a desired position
TODO: use magnetic field values and shower core from config-file
'''
import numpy as np
from scipy import signal
import operator
import logging
import os

from os.path import split
import sys
from copy import deepcopy
from frame import UVWGetter

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

import hdf5fileinout as hdf5io

import astropy.units as u


#======================================
def unwrap(phi, ontrue=None):
    """Unwrap the phase so that the absolute difference
      between 2 consecutive phases remains below Pi

    Parameters:
    ----------
        phi: numpy array, float
            phase of the signal trace
        ontrue: str
            printing option, default=None

    Returns:
    ----------
        phi_unwrapped: numpy array, float
            unwarpped phase of the signal trace

    Adapted by E. Hivon (2020-02) from A. Zilles' unwrap
    """
    eps = np.finfo(np.pi).resolution
    thr = np.pi - eps
    pi2 = 2. * np.pi
    phi_unwrapped = np.zeros(phi.shape)
    p0  = phi_unwrapped[0] = phi[0]
    l   = 0
    for i0, p1 in enumerate(phi[1:]):
        i = i0 + 1
        dp = p1 - p0
        if (np.abs(dp) > thr):
            dl = np.floor_divide(abs(dp), pi2) + 1
            if (dp > 0):
                l -= dl
            else:
                l += dl
        phi_unwrapped[i] = p1 + l * pi2
        p0 = p1
        if ontrue is not None:
            print(i, phi[i],           phi[i-1],           abs(phi[i] - phi[i-1]),
                  l, phi_unwrapped[i], phi_unwrapped[i-1], abs(phi_unwrapped[i] - phi_unwrapped[i-1]))

    return phi_unwrapped
#======================================


def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):

	ndotu = planeNormal.dot(rayDirection)
	if abs(ndotu) < epsilon:
		raise RuntimeError("no intersection or line is within plane")

	w = rayPoint - planePoint
	si = -planeNormal.dot(w) / ndotu
	Psi = w + si * rayDirection + planePoint
	return Psi


    ##Define plane
	#planeNormal = np.array([0, 0, 1])
	#planePoint = np.array([0, 0, 5]) #Any point on the plane

	##Define ray
	#rayDirection = np.array([0, -1, -1])
	#rayPoint = np.array([0, 0, 10]) #Any point along the ray

	#Psi = LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint)
	#print ("intersection at", Psi)





def interpolate_trace(t1, trace1, x1, t2, trace2, x2, xdes, upsampling=None,  zeroadding=None):
    """Interpolation of signal traces at the specific position in the frequency domain

    The interpolation of traces needs as input antenna position 1 and 2, their traces (filtered or not)
    in one component, their time, and the desired antenna position and returns the trace ( in x,y,z coordinate system) and the time from the desired antenna position.
    Zeroadding and upsampling of the signal are optional functions.

    IMPORTANT NOTE:
    The interpolation of the phases includes the interpolation of the signal arrival time. A linear interpolation implies a plane radio
    emission wave front, which is a simplification as it is hyperbolic in shape. However, the wave front can be estimated as a plane between two simulated observer positions
    for a sufficiently dense grid of observers, as then parts of the wave front are linear on small scales.

    This script bases on the diploma thesis of Ewa Holt (KIT, 2013) in the context of AERA/AUGER. It is based on the interpolation of the amplitude and the pahse in the frequency domain.
    This can lead to misidentifiying of the correct phase. We are working on the interplementaion on a more robust interpolation of the signal time.
    Feel free to include it if you have some time to work on it. The script is completely modular so that single parts can be substitute easily.


    Parameters:
    ----------
            t1: numpy array, float
                time in ns of antenna 1
            trace1: numpy array, float
                single component of the electric field's amplitude of antenna 1
            x1: numpy array, float
                position of antenna 1
            t2: numpy array, float
                time in ns of antenna 2
            trace2: numpy array, float
                single component of the electric field's amplitude of antenna 2
            x2: numpy array, float
                position of antenna 2
            xdes: numpy arry, float
                antenna position for which trace is desired, in meters
            upsampling: str
                optional, True/False, performs upsampling of the signal, by a factor 8
            zeroadding: str
                optional, True/False, adds zeros at the end of the trace of needed

    Returns:
    ----------
        xnew: numpy array, float
            rough estimation of the time for signal at desired antenna position in ns
        tracedes: numpy array, float
            interpolated electric field component at desired antenna position
    """
    DISPLAY = False
    ontrue = None

    # hand over time traces of one efield component -t1=time, trace1=efield- and the position
    # x1 of the first antenna, the same for the second antenna t2,trace2, x2.
    # xdes is the desired antenna position (m) where you would like to have the efield trace in time
    # if necessary you have to do an upsampling of the trace: upsampling=On

    factor_upsampling = 1
    if upsampling is not None:
        factor_upsampling = 8

    # calculating weights: should be done with the xyz coordinates
    # since in star shape pattern it is mor a radial function connection the poistion of
    # same signal as linear go for that solution.
    # if lines ar on a line, it will give the same result as before
    tmp1 = np.linalg.norm(x1 - xdes)
    tmp2 = np.linalg.norm(x2 - xdes)

    tmp = 1. / (tmp1 + tmp2)
    weight1 = tmp2 * tmp
    weight2 = tmp1 * tmp

    if np.isinf(weight1):
        print("weight = inf")
        print(x1, x2, xdes)
        weight1 = 1.
        weight2 = 0.
    if np.isnan(weight1):
        print('Attention: projected positions equivalent')
        weight1 = 1.
        weight2 = 0.
    epsilon = np.finfo(float).eps
    if (weight1 > 1. + epsilon) or (weight2 > 1 + epsilon):
        print("weight larger 1: ", weight1, weight2, x1, x2, xdes, np.linalg.norm(
            x2-x1), np.linalg.norm(x2-xdes), np.linalg.norm(xdes-x1))
    if weight1 + weight2 > 1 + epsilon:
        print("PulseShape_Interpolation.py: order in simulated positions. Check whether ring or ray structure formed first")
        print(weight1, weight2, weight1 + weight2)


    #################################################################################
    # Fourier Transforms

    # first antenna
    # upsampling if necessary
    if upsampling is not None:
        trace1 = signal.resample(trace1, len(trace1)*factor_upsampling)
        t1 = np.linspace(t1[0], t1[-1], len(trace1)
                            * factor_upsampling, endpoint=False)

    if zeroadding is True:
        max_element = len(trace1)  # to shorten the array after zeroadding
        print(max_element)
        xnew = np.linspace(t1[0], 1.01*t1[-1],
                              int((1.01*t1[-1]-t1[0])/(t1[2]-t1[1])))
        print(len(xnew))
        #xnew = xnew*1.e-9  # ns -> s
        zeros = np.zeros(len(xnew)-max_element)
        f = trace1
        f = np.hstack([f, zeros])

    if zeroadding is None:
        f = trace1
        xnew = t1 #*1.e-9

    fsample = 1./((xnew[1]-xnew[0]))  # GHz

    freq = np.fft.rfftfreq(len(xnew), 1./fsample)
    FFT_Ey = np.fft.rfft(f)

    Amp = np.abs(FFT_Ey)
    phi = np.angle(FFT_Ey)
    #Eric
    phi_unwrapped = unwrap(phi, ontrue)

    #############################

    # second antenna
    # t in ns, Ex in muV/m, Ey, Ez
    # NOTE: Time binning always 1ns

    # upsampling if needed
    if upsampling is not None:
        trace = signal.resample(trace2, len(trace2)*factor_upsampling)
        trace2 = trace
        t2 = np.linspace(t2[0], t2[-1], len(trace2)
                            * factor_upsampling, endpoint=False)

    if zeroadding is True:
        # get the same length as xnew
        xnew2 = np.linspace(
            t2[0], t2[0] + (xnew[-1]-xnew[0])*1e9, len(xnew))
        #xnew2 = xnew2*1.e-9
        f2 = trace2
        f2 = np.hstack([f2, zeros])

    if zeroadding is None:
        f2 = trace2
        xnew2 = t2 #*1e-9  # ns -> s

    fsample2 = 1./((xnew2[1]-xnew2[0]))  # *1.e-9 to get time in s

    freq2 = np.fft.rfftfreq(len(xnew2), 1./fsample2)
    FFT_Ey = np.fft.rfft(f2)

    Amp2 = np.abs(FFT_Ey)
    phi2 = np.angle(FFT_Ey)
    #Eric
    phi2_unwrapped = unwrap(phi2, ontrue)

    ### Get the pulse sahpe at the desired antenna position

    # get the phase
    #Eric
    phides = weight1 * phi_unwrapped + weight2 * phi2_unwrapped

    if ontrue is not None:
        print(phides)
    if DISPLAY:
        phides2 = phides.copy()

    #Eric re-unwrap: get -pi to +pi range back and check whether phidesis in between (im not wraping any more)
    phides = np.mod(phides + np.pi, 2. * np.pi) - np.pi

    #################################################################################
    ### linearly interpolation of the amplitude

    #Amp, Amp2
    # Since the amplitude shows a continuous unipolar shape, a linear interpolation is sufficient

    Ampdes = weight1 * Amp + weight2 * Amp2
    if DISPLAY:
        Ampdes2 = Ampdes.copy()

    # inverse FFT for the signal at the desired position
    Ampdes = Ampdes.astype(np.complex64)
    phides = phides.astype(np.complex64)
    if DISPLAY:
        phides2 = phides2.astype(np.complex64)
    Ampdes *= np.exp(1j * phides)

    tracedes = (np.fft.irfft(Ampdes))
    tracedes = tracedes.astype(float)

    #this is a crude interpolation of the time
    tdes=(xnew*weight1+xnew2*weight2)
    if(len(tdes)>len(tracedes)):
     tdes=tdes[0:-1] #and this is required becouse the inverse fft returns one less time bin
    #print("weights 1:"+str(weight1)+ " 2:"+str(weight2))


    # PLOTTING

    if (DISPLAY):
        import matplotlib.pyplot as plt
        import pylab

        fig1 = plt.figure(1, dpi=120, facecolor='w', edgecolor='k')
        plt.subplot(211)
        plt.plot(freq, phi, 'ro-', label="first")
        plt.plot(freq2, phi2, 'bo-', label="second")
        plt.plot(freq2, phides, 'go--', label="interpolated")
        #plt.plot(freq2, phi_test, 'co--', label= "real")
        plt.xlabel(r"Frequency (GHz)", fontsize=16)
        plt.ylabel(r"phase (rad)", fontsize=16)
        #plt.xlim(flow, fhigh)
        print(len(freq), len(phi),len(Amp), len(freq2),len(phi2),len(Amp2),len(phides),len(Ampdes2))

        #pylab.legend(loc='upper left')

        #plt.subplot(312)
        #ax = fig1.add_subplot(3, 1, 2)
        #plt.plot(freq, phi_unwrapped, 'r+')
        #plt.plot(freq2, phi2_unwrapped, 'bx')
        #plt.plot(freq2, phides2, 'g^')
        #plt.plot(freq2, phi_test_unwrapped, 'c^')
        #plt.xlabel(r"Frequency (Hz)", fontsize=16)
        #plt.ylabel(r"phase (rad)", fontsize=16)
        # plt.show()
        # plt.xlim([0,0.1e8])
        # plt.xlim([1e8,2e8])
        # plt.ylim([-10,10])
        # ax.set_xscale('log')
        #plt.xlim(flow, fhigh)

        plt.subplot(212)
        plt.plot(freq, Amp, 'r+')
        plt.plot(freq2, Amp2, 'bx')
        plt.plot(freq2, Ampdes2, 'g^')
        #plt.plot(freq2, Amp_test, 'c^')
        plt.xlabel(r"Frequency (GHz)", fontsize=16)
        plt.ylabel(r"Amplitude muV/m/GHz ", fontsize=16)
        #print("Min Amplitude: " + str(np.min(Amp)) + " Amplitude 2: " + str(np.min(Amp2)))
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        #plt.ylim([1e1, 10e3])
        #plt.xlim(flow, fhigh)

        plt.show()

################################## CONTROL

    if DISPLAY:
        ##### PLOTTING

        import matplotlib.pyplot as plt
        plt.plot(np.real(t1), np.real(trace1), 'g:', label= "antenna 1")
        plt.plot(np.real(t2), np.real(trace2), 'b:', label= "antenna 2")
        plt.plot(np.real(tdes), np.real(tracedes), 'r-', label= "Synthetized")

        plt.xlabel(r"time (ns)", fontsize=16)
        plt.ylabel(r"Amplitude muV/m ", fontsize=16)
        plt.legend(loc='best')

        plt.show()

    if zeroadding is True:
        # hand over time of first antenna since interpolation refers to that time
        return tdes[0:max_element], tracedes[0:max_element]

    if upsampling is not None:
        return tdes[0:-1:8], tracedes[0:-1:8]
    else:
        #xnew = np.delete(xnew, -1)
        return tdes, tracedes  # back to ns



def do_interpolation_hdf5(desired, InputFilename, OutputFilename, antennamin=0, antennamax=159, EventNumber=0, DISPLAY=False, usetrace='efield'):
    '''
    Reads in arrays, looks for neighbours, calls the interpolation and saves the traces

    Parameters:
    ----------
    desired: str
        numpy array of desired antenna positions (x,y,z,t0 info)
    InputFilename: str
        path to HDF5 simulation file
        The script accepts starshape as well as grid arrays
    antennamin,antennamax:int
        the program is designed to run on the first 160 antennas. If your simulation has more, you can specify a range to be used...but it has to be tested
    EventNumber: int
        number of event in the file to use. you can process only one at a time
    shower_core: numpy array
        position of shower core for correction (NOT TESTED, CHECK BEFORE USING IT!) In ZHAireS, its the altitude of the simulation, but i dont really get why.
    DISPLAY: True/False
        enables printouts and plots
    usetrace: str (note that for now you can only do one at a time, and on different output files)
        efield
        voltage
        filteredvoltage

    Returns:
    ----------
        --
    Saves traces via index infomation in same folder as desired antenna positions


    NOTE: The selection of the neigbours is sufficiently stable, but does not always pick the "best" neigbour, still looking for an idea
    '''
    #print(shower_core)
    DEVELOPMENT=True #only for developing, use when working on the starshape patern and trying to interpolate the random check antenas on the starshape (or it will crash)
                     #it disables removing antennas outside the pattern

    #0 dont plot
    #3 plot starshape with selected antennas and interpolated signals, and interpolated signals errors if in DEVELOPMENT
    #4 plot antenna positions in ground and in UVW
    #5 plot antenna cuadrants and selection

    DISPLAY=0

    projection="Conical"

    PLOTPAPER=True

    print("using projection ", projection)
    starshape_exploit=True

    if(starshape_exploit): #If i found the two, "external" points, then the internals must be the same minus 8. This is only valid in the 8 branches starshape we use
      print("exploiting starshape for interpolation, this is experimental...use with extreme criticism!")


    CurrentEventNumber=EventNumber
    CurrentRunInfo=hdf5io.GetRunInfo(InputFilename)
    CurrentNumberOfEvents=hdf5io.GetNumberOfEvents(CurrentRunInfo)
    CurrentEventName=hdf5io.GetEventName(CurrentRunInfo,CurrentEventNumber)
    CurrentEventInfo=hdf5io.GetEventInfo(InputFilename,CurrentEventName)
    CurrentShowerSimInfo=hdf5io.GetShowerSimInfo(InputFilename,CurrentEventName)
    CurrentSignalSimInfo=hdf5io.GetSignalSimInfo(InputFilename,CurrentEventName)

    Zenith=hdf5io.GetEventZenith(CurrentRunInfo,0)
    Azimuth=hdf5io.GetEventAzimuth(CurrentRunInfo,0)


    PhiGeo= hdf5io.GetEventBFieldDecl(CurrentEventInfo)
    PhiGeo= 0 #we are asuming magnetic coordinates
    print("we are assuming azimuths are magnetic!")
    ThetaGeo= hdf5io.GetEventBFieldIncl(CurrentEventInfo) + 90.0 #adjust to GRAND coordinates.
    GroundAltitude=hdf5io.GetGroundAltitude(CurrentEventInfo)
    shower_core=np.array([0,0,0])#np.array([0,0,GroundAltitude])
    XmaxPosition=hdf5io.GetXmaxPosition(CurrentEventInfo)
    XmaxDistance=hdf5io.GetEventXmaxDistance(CurrentRunInfo,0)

    tbinsize=hdf5io.GetTimeBinSize(CurrentSignalSimInfo)
    tmin=hdf5io.GetTimeWindowMin(CurrentSignalSimInfo)
    tmax=hdf5io.GetTimeWindowMax(CurrentSignalSimInfo)

    # SIMULATION
    # Read in simulated position list
    CurrentAntennaInfo=hdf5io.GetAntennaInfo(InputFilename,CurrentEventName)
    antennamax=antennamax+1
    #one way of putting the antenna information as the numpy array this script was designed to use:
    xpoints=CurrentAntennaInfo['X'].data[antennamin:antennamax]
    ypoints=CurrentAntennaInfo['Y'].data[antennamin:antennamax]
    zpoints=CurrentAntennaInfo['Z'].data[antennamin:antennamax]
    t0s=CurrentAntennaInfo['T0'].data[antennamin:antennamax]
    positions_sims=np.column_stack((xpoints,ypoints,zpoints))

    # DESIRED
    # Hand over a list file including the antenna positions you would like to have. This could be improved by including an ID.
    positions_des = desired #np.loadtxt(desired,usecols=(2,3,4))


    if DISPLAY>0:
        print('desired positions: '+ str(len(positions_des)))
        #print(positions_des, len(positions_des))
    if len(positions_des) <=1:
        print("Files of desired positions has to consist of at least two positions, Bug to be fixed")

    if DISPLAY>0:
        print('simulated positions: ' + str(len(positions_sims)))
        #print(positions_sims, len(positions_sims))
    if len(positions_sims) <=1:
        print("Files of simulated positions has to consist of at least two positions, Bug to be fixed")

    #making the table of desired antennas for the file (but we will save it later, when we know wich were actually used)
    DesiredAntennaInfoMeta=hdf5io.CreatAntennaInfoMeta(split(InputFilename)[1],CurrentEventName,AntennaModel="Interpolated")
    DesiredIds=np.arange(0, len(positions_des)) #this could be taken from the input file of desired antennas
    DesiredAntx=deepcopy(positions_des.T[0])
    DesiredAnty=deepcopy(positions_des.T[1])
    DesiredAntz=deepcopy(positions_des.T[2]) #this deepcopy bullshit is becouse position_des is later modified by the rotation, and transposition apparently creates a shallow copy (a reference)
    DesiredSlopeA=np.zeros(len(positions_des))
    DesiredSlopeB=np.zeros(len(positions_des))
    DesiredT0=deepcopy(positions_des.T[3])

    #now i come back to having only x,y,z on positions_des
    positions_des=desired[:,0:3]

    #------------------------Plot in XYZ Coord
    if DISPLAY>3:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        #ax = fig.gca(projection='3d')

        #ax.scatter(positions_sims[:,0], positions_sims[:,1], positions_sims[:,2], label = "Simulated")
        #ax.scatter(positions[:,0], positions[:,1], positions[:,2], label = "Synthetized")
        #ax.scatter(shower_core[0], shower_core[1], shower_core[2], label = "Shower core")
        for j in range(0,len(positions_des[:,1])):
                ax.annotate(str(j), ((positions_des[j,0], positions_des[j,1])))
        ax.scatter(positions_sims[:,0], positions_sims[:,1], label = "Simulated")
        ax.scatter(positions_des[:,0], positions_des[:,1], label = "Synthetized")
        ax.scatter(shower_core[0], shower_core[1],  label = "Shower core")

        plt.title("XYZ coordinates")
        plt.legend(loc=2)
        plt.show()
    #------------------------


    ##--##--##--##--##--##--##--##--##-##--##--##-##--##--## START: WRAP UP AS FUNCTION (PROJECTION AND ROTATION)
    #### START: UNDO projection
    #define shower vector
    az_rad=np.deg2rad(180.+Azimuth)#Note ZHAIRES units used
    zen_rad=np.deg2rad(180.-Zenith)

    # shower vector  = direction of line for backprojection, TODO should be substituded bey line of sight Xmax - positions
    v = np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])
    v = v/np.linalg.norm(v)

    pos_des= np.zeros([len(positions_des[:,1]),3])
    pos_sims= np.zeros([len(positions_sims[:,1]),3])
    #GetUVW = UVWGetter(shower_core[0], shower_core[1], shower_core[2], zenith, azimuth, phigeo, thetageo)
    GetUVW = UVWGetter(0., 0., 0., Zenith, Azimuth, PhiGeo, ThetaGeo)

    if(projection=="Geometric"):
      ### START: ROTATE INTO SHOWER COORDINATES, and core for offset by core position, alreadz corrected in projection
      # for back projection position vector line is projected position
      # for back projection normal vector of plane to intercsect == v
      n = v

      for i in np.arange(0,len(positions_des[:,1])):
        b=-np.dot(n,positions_des[i,:])/ np.dot(n, v)
        positions_des[i,:] = positions_des[i,:] + b*v - shower_core # correct by shower core position
      for i in np.arange(0,len(positions_sims[:,1])):
        b=-np.dot(n,positions_sims[i,:])/ np.dot(n, v)
        positions_sims[i,:] = positions_sims[i,:] + b*v - shower_core # correct by shower core position

      #rotate desired positions
      for i in np.arange(0,len(positions_des[:,1])):
         pos_des[i,:]=GetUVW(positions_des[i,:], )

      # Rotate simulated positions
      for i in np.arange(0,len(positions_sims[:,1])):
        pos_sims[i,:] = GetUVW(positions_sims[i,:], )

    if(projection=="Conical"):

      planeNormal=v
      XmaxPosition=v*(XmaxDistance+3000) #origin on 0,0,ground
      print("asuming cone vertex is 3Km behind Xmax")
      rayPoint=XmaxPosition
      planePoint=np.array([0,0,0]) #the starshape is always on the ground when generated for ZHAireS

      #rotate desired positions
      for i in np.arange(0,len(positions_des[:,1])):
         #print("Antenna",i)
         rayDirection=positions_des[i,:]-np.array([0,0,GroundAltitude])-XmaxPosition
         #print("rayDirection",rayDirection)
         pos_des[i,:]=LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6)
         #print("collision",pos_des[i,:])
         pos_des[i,:]=GetUVW(pos_des[i])# as used later to go in vxB
         #print("UVW",pos_des[i,:])

      #rotate simulation positions
      for i in np.arange(0,len(positions_sims[:,1])):
         #print("Antenna",i)
         rayDirection=positions_sims[i,:]-np.array([0,0,GroundAltitude])-XmaxPosition
         #print("rayDirection",rayDirection)
         pos_sims[i,:]=LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6)
         #print("collision",pos_sims[i,:])
         pos_sims[i,:]=GetUVW(pos_sims[i])# as used later to go in vxB
         #print("UVW",pos_sims[i,:])

    # ------------------ PLot in UVW coordinates
    if DISPLAY>3:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(1,1,1)
        #ax2 = fig2.gca(projection='3d')
        #ax2.scatter(pos_sims[:,0], pos_sims[:,1], pos_sims[:,2], label = "Simulated")
        #ax2.scatter(pos_des[:,0], pos_des[:,1], pos_des[:,2], label = "Synthetized")
        for j in range(0,len(pos_des[:,1])):
                ax2.annotate(str(j), ((pos_des[j,1], pos_des[j,2])))

        ## x component should be 0
        ax2.scatter(pos_sims[:,1], pos_sims[:,2], label = "Simulated")
        ax2.scatter(pos_des[:,1], pos_des[:,2], label = "Synthetized")
        ax2.scatter(0, 0, marker ="x",  label = "core")

        plt.title("shower coordinates")
        plt.legend(loc=2)
        plt.show()
    # ------------------

    # calculate radius and angle for simulated positions and store some in list
    points=[]

    for i in np.arange(0,len(pos_sims[:,1])):  # position should be within one plane yz plane, remove x=v component for simplicity
        #points.append([i, pos_sims[i,1], pos_sims[i,2] ])
        theta2 = np.arctan2(pos_sims[i,2], pos_sims[i,1])
        radius2 = np.sqrt( pos_sims[i,1]**2 + pos_sims[i,2]**2 )
        if round(theta2,4) == -3.1416:        #this is to put theta 2 in (-3.1416,3.1416)
            theta2*=-1
        points.append([i, theta2, radius2])

    # Here the magic starts, so i will write the file header

    #not using them, but i put SignalSim And ShowerSim Info
    #For now, i save a copy. I could modify some fields to show this is an interpolation
    hdf5io.SaveRunInfo(OutputFilename,CurrentRunInfo)
    hdf5io.SaveEventInfo(OutputFilename,CurrentEventInfo,CurrentEventName)
    hdf5io.SaveShowerSimInfo(OutputFilename,CurrentShowerSimInfo,CurrentEventName)
    hdf5io.SaveSignalSimInfo(OutputFilename,CurrentSignalSimInfo,CurrentEventName)

    #this is ugly, becouse i repeat the determination of the closest antenna. This could be done only once...but my time is limited too!
    if(usetrace=="all"):
      print("usetrace is all, looping over all trace types")
      usetracelist=["efield","voltage","filteredvoltage"]
    else:
      usetracelist=[str(usetrace)]

    for tracetype in usetracelist:
        print("computing for "+tracetype)
        remove_antenna=[] #this will be cleared for each trace type we loop on, but since the 3 are discarding the same antenas, its ok

        #i loops only over desired in-plane positions, acting as new reference
        for i in np.arange(0,len(pos_des[:,1])):  # position should be within one plane yz plane, remove x=v component for simplicity
            #print("desired antena:"+str(i))
            theta = np.arctan2(pos_des[i,2], pos_des[i,1])
            radius = np.sqrt( pos_des[i,1]**2 + pos_des[i,2]**2 )
            #print("index of desired antenna ", ind[i], theta, radius, )


            # The 4 quadrants -- in allen 4 Ecken soll Liebe drin stecken
            points_I=[]
            points_II=[]
            points_III=[]
            points_IV=[]

            #m loops over the simulated positions
            for m in np.arange(0,len(points)): # desired points as reference
                delta_phi = points[m][1]-theta
                if delta_phi > np.pi:
                    delta_phi = delta_phi -2.*np.pi
                elif delta_phi < -np.pi:
                    delta_phi = delta_phi + 2.*np.pi


                delta_r = points[m][2]-radius

                #distance = np.sqrt(delta_r**2 + (delta_r *delta_phi)**2 ) # weighting approach1
                #distance= np.sqrt((pos_sims[m,1]-pos_des[i,1])**2. +(pos_sims[m,2]-pos_des[i,2])**2.) # euclidean distance
                distance= np.sqrt(points[m][2]**2. +radius**2. -2.*points[m][2]*radius* np.cos(points[m][1]-theta) ) #polar coordinates

                if delta_phi >= 0. and  delta_r >= 0:
                    points_I.append((m,delta_phi,delta_r, distance))
                if delta_phi >= 0. and  delta_r <= 0:
                    points_II.append((m,delta_phi,delta_r, distance))
                if delta_phi <= 0. and  delta_r <= 0:
                    points_III.append((m,delta_phi,delta_r, distance))
                if delta_phi <= 0. and  delta_r >= 0:
                    points_IV.append((m,delta_phi,delta_r, distance))

            bailoutI=0
            bailoutII=0
            bailoutIII=0
            bailoutIV=0
            if not points_I:
                #print("list - Quadrant 1 - empty --> no interpolation for ant", str(ind[i]))
                bailoutI=1
            if not points_II:
                #print("list - Quadrant 2 - empty --> no interpolation for ant", str(ind[i]))
                bailoutII=1
                points_II=points_I
            if not points_III:
                #print("list - Quadrant 3 - empty --> no interpolation for ant", str(ind[i]))
                bailoutIII=1
                points_III=points_IV
            if not points_IV:
                #print("list - Quadrant 4 - empty --> no interpolation for ant", str(ind[i]))
                bailoutIV=1

            if(bailoutII==1 and bailoutIII==1 and bailoutIV==0 and bailoutI==0):
                print("I cannot find antennas with lower r, lets try using only the inner antennas")

            if(bailoutI==1 or bailoutIV==1 or (bailoutII==1 and bailoutIII==0) or (bailoutII==0 and bailoutII==1)):
              #antenas to remove
              print(str(i)+" antenna is outside of the starshape patern, removing it",pos_des[i,:],theta,radius)
              if(DEVELOPMENT==False):
                remove_antenna.append(i)
              else:
               # all this crap for putting empty traces is becouse im not ready to process files that only have some traces.

               AntennaID=hdf5io.GetAntennaID(CurrentAntennaInfo,0) #this is just to get the sape of what i have to save with 0s
               if(tracetype=='efield'):
                  txt0=hdf5io.GetAntennaEfield(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
               elif(tracetype=='voltage'):
                  txt0=hdf5io.GetAntennaVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
               elif(tracetype=='filteredvoltage'):
                  txt0=hdf5io.GetAntennaFilteredVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
               else:
                   print("You must specify either efield, voltage or filteredvoltage, bailing out")

               if(tracetype=='efield'):
                 efield=np.zeros(np.shape(txt0))
                 EfieldTable=hdf5io.CreateEfieldTable(efield, CurrentEventName, CurrentEventNumber , DesiredIds[i], i, "Interpolated", info={})
                 hdf5io.SaveEfieldTable(OutputFilename,CurrentEventName,str(DesiredIds[i]),EfieldTable)
               elif(tracetype=='voltage'):
                 voltage=np.zeros(np.shape(txt0))
                 VoltageTable=hdf5io.CreateVoltageTable(voltage, CurrentEventName, CurrentEventNumber , DesiredIds[i], i, "Interpolated", info={})
                 hdf5io.SaveVoltageTable(OutputFilename,CurrentEventName,str(DesiredIds[i]),VoltageTable)
               elif(tracetype=='filteredvoltage'):
                 filteredvoltage=np.zeros(np.shape(txt0))
                 VoltageTable=hdf5io.CreateVoltageTable(filteredvoltage, CurrentEventName, CurrentEventNumber , DesiredIds[i], i, "Interpolated", info={})
                 hdf5io.SaveFilteredVoltageTable(OutputFilename,CurrentEventName,str(DesiredIds[i]),VoltageTable)

            else:
                # ------------------

                points_I=np.array(points_I, dtype = [('index', 'i4'), ('delta_phi', 'f4'), ('delta_r', 'f4'), ('distance', 'f4')])
                points_II=np.array(points_II, dtype = [('index', 'i4'), ('delta_phi', 'f4'), ('delta_r', 'f4'), ('distance', 'f4')])
                points_III=np.array(points_III, dtype = [('index', 'i4'), ('delta_phi', 'f4'), ('delta_r', 'f4'), ('distance', 'f4')])
                points_IV=np.array(points_IV, dtype = [('index', 'i4'), ('delta_phi', 'f4'), ('delta_r', 'f4'), ('distance', 'f4')])

                ## Sort points; not optimal (the best) solution for all, but brings stable/acceptable results
                points_I = np.sort(points_I, order=['distance', 'delta_phi', 'delta_r'])
                points_II = np.sort(points_II, order=['distance', 'delta_phi', 'delta_r'])
                points_III = np.sort(points_III, order=[ 'distance','delta_phi', 'delta_r'])
                points_IV = np.sort(points_IV, order=['distance', 'delta_phi', 'delta_r'])
                #indices of 4 closest neigbours: points_I[0][0], points_II[0][0], points_III[0][0], points_IV[0][0]

                # try to combine the one with roughly the same radius first and then the ones in phi
                # since the weighting is done just for the distance, and that is the only thing why the position is used inside that function,
                #i will input to interpolate_trace the y and z components, and make x=0 (this is becouse in the uxVxB plane one component is 0?

                if(starshape_exploit): #If i found the two, "external" points, then the internals must be the same minus 8. This is only valid in the 8 branches starshape we use
                  if(points_I[0][0]>7):
                    record=np.where(points_II['index']==points_I[0][0]-8)
                    points_II[0]=points_II[record]
                  if(points_IV[0][0]>7):
                    record=np.where(points_III['index']==points_IV[0][0]-8)
                    points_III[0]=points_III[record]

                #points I and IV have a higher radius. lets take the average radius, and the same theta of the desired point as the point1
                #this is stored alreadty in the opoints variable
                meanr=(points[points_I[0][0]][2]+points[points_IV[0][0]][2])/2.0
                #and theta is already available
                point_online1=np.array([0,meanr*np.cos(theta),meanr*np.sin(theta)])

                meanr=(points[points_II[0][0]][2]+points[points_III[0][0]][2])/2.0
                #and theta is already available
                point_online2=np.array([0,meanr*np.cos(theta),meanr*np.sin(theta)])
                #points II and III have a lower radius. lets take the average radius, and the same theta of the desired point as the point2

                #this is to show a list of the indices of antennas in the each quadrant
                if (DISPLAY>4):
                  listI=list(list(zip(*points_I))[0])
                  listII=list(list(zip(*points_II))[0])
                  listIII=list(list(zip(*points_III))[0])
                  listIV=list(list(zip(*points_IV))[0])

                  mypoints_I=[]
                  mypoints_II=[]
                  mypoints_III=[]
                  mypoints_IV=[]

                  for h in listI:
                    mypoints_I.append((pos_sims[h,1],pos_sims[h,2]))

                  for h in listII:
                    mypoints_II.append((pos_sims[h,1],pos_sims[h,2]))

                  for h in listIII:
                    mypoints_III.append((pos_sims[h,1],pos_sims[h,2]))

                  for h in listIV:
                    mypoints_IV.append((pos_sims[h,1],pos_sims[h,2]))

                  mypoints_I=np.array(mypoints_I)
                  mypoints_II=np.array(mypoints_II)
                  mypoints_III=np.array(mypoints_III)
                  mypoints_IV=np.array(mypoints_IV)

                  fig3a = plt.figure()
                  ax3a = fig3a.add_subplot(1,1,1)

                  for j in range(0,len(pos_sims[:,1])):
                    ax3a.annotate(str(j), ((pos_sims[j,1], pos_sims[j,2])))

                  ## x component should be 0
                  ax3a.scatter(pos_des[i,1], pos_des[i,2], label = "Synthetized")
                  ax3a.scatter(mypoints_I[:,0], mypoints_I[:,1], s=90,marker ="D",label = "1")
                  ax3a.scatter(mypoints_II[:,0], mypoints_II[:,1], s=90, marker ="D",label = "2")
                  ax3a.scatter(mypoints_III[:,0], mypoints_III[:,1], s=90, marker ="D",label = "3")
                  ax3a.scatter(mypoints_IV[:,0], mypoints_IV[:,1], s=90,marker ="D",label = "4")

                  ax3a.scatter(pos_sims[points_I[0][0],1], pos_sims[points_I[0][0],2], s=20, marker ="D",label = "1")
                  ax3a.scatter(pos_sims[points_II[0][0],1], pos_sims[points_II[0][0],2], s=20, marker ="D", label = "2")
                  ax3a.scatter(pos_sims[points_III[0][0],1], pos_sims[points_III[0][0],2], s=20, marker ="D", label = "3")
                  ax3a.scatter(pos_sims[points_IV[0][0],1], pos_sims[points_IV[0][0],2], s=20, marker ="D",label = "4")

                  ax3a.scatter(point_online1[1], point_online1[2], marker ="x")
                  ax3a.scatter(point_online2[1], point_online2[2], marker ="x")
                  ax3a.scatter(0, 0, marker ="x",  label = "core")
                  ax3a.plot([0, pos_des[i,1]], [0, pos_des[i,2]])
                  plt.title("my"+str(i))
                  plt.legend(loc=2)
                  plt.show()
                # ------------------

                ## the interpolation of the pulse shape is performed, in x, y and z component
                AntennaID=hdf5io.GetAntennaID(CurrentAntennaInfo,points_I[0][0])
                if(tracetype=='efield'):
                  txt0=hdf5io.GetAntennaEfield(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='voltage'):
                  txt0=hdf5io.GetAntennaVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='filteredvoltage'):
                  txt0=hdf5io.GetAntennaFilteredVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                else:
                  print('You must specify either efield, voltage or filteredvoltage, bailing out')
                  return 0

                AntennaID=hdf5io.GetAntennaID(CurrentAntennaInfo,points_IV[0][0])
                if(tracetype=='efield'):
                  txt1=hdf5io.GetAntennaEfield(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='voltage'):
                  txt1=hdf5io.GetAntennaVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='filteredvoltage'):
                  txt1=hdf5io.GetAntennaFilteredVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                else:
                  print('You must specify either efield, voltage or filteredvoltage, bailing out')
                  return 0

                AntennaID=hdf5io.GetAntennaID(CurrentAntennaInfo,points_II[0][0])
                if(tracetype=='efield'):
                  txt2=hdf5io.GetAntennaEfield(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='voltage'):
                  txt2=hdf5io.GetAntennaVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='filteredvoltage'):
                  txt2=hdf5io.GetAntennaFilteredVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                else:
                  print('You must specify either efield, voltage or filteredvoltage, bailing out')
                  return 0


                AntennaID=hdf5io.GetAntennaID(CurrentAntennaInfo,points_III[0][0])
                if(tracetype=='efield'):
                  txt3=hdf5io.GetAntennaEfield(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='voltage'):
                  txt3=hdf5io.GetAntennaVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                elif(tracetype=='filteredvoltage'):
                  txt3=hdf5io.GetAntennaFilteredVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                else:
                  print('You must specify either efield, voltage or filteredvoltage, bailing out')
                  return 0

                matias_patch_to_point_I=np.array([0,pos_sims[points_I[0][0]][1],pos_sims[points_I[0][0]][2]])
                matias_patch_to_point_II=np.array([0,pos_sims[points_II[0][0]][1],pos_sims[points_II[0][0]][2]])
                matias_patch_to_point_III=np.array([0,pos_sims[points_III[0][0]][1],pos_sims[points_III[0][0]][2]])
                matias_patch_to_point_IV=np.array([0,pos_sims[points_IV[0][0]][1],pos_sims[points_IV[0][0]][2]])

                xnew1, tracedes1x = interpolate_trace(txt0.T[0], txt0.T[1], matias_patch_to_point_I , txt1.T[0], txt1.T[1], matias_patch_to_point_IV, point_online1 ,upsampling=None, zeroadding=None)
                xnew1, tracedes1y = interpolate_trace(txt0.T[0], txt0.T[2], matias_patch_to_point_I , txt1.T[0], txt1.T[2], matias_patch_to_point_IV, point_online1 ,upsampling=None, zeroadding=None)
                xnew1, tracedes1z = interpolate_trace(txt0.T[0], txt0.T[3], matias_patch_to_point_I , txt1.T[0], txt1.T[3], matias_patch_to_point_IV, point_online1 ,upsampling=None, zeroadding=None)

                xnew2, tracedes2x = interpolate_trace(txt2.T[0], txt2.T[1], matias_patch_to_point_II , txt3.T[0], txt3.T[1], matias_patch_to_point_III, point_online2 ,upsampling=None, zeroadding=None)
                xnew2, tracedes2y = interpolate_trace(txt2.T[0], txt2.T[2], matias_patch_to_point_II , txt3.T[0], txt3.T[2], matias_patch_to_point_III, point_online2 ,upsampling=None, zeroadding=None)
                xnew2, tracedes2z = interpolate_trace(txt2.T[0], txt2.T[3], matias_patch_to_point_II , txt3.T[0], txt3.T[3], matias_patch_to_point_III, point_online2 ,upsampling=None, zeroadding=None)

                ###### Get the pulse shape of the desired position from projection on line1 and 2
                xnew_desiredx, tracedes_desiredx =interpolate_trace(xnew1, tracedes1x, point_online1, xnew2, tracedes2x, point_online2, np.array([0,pos_des[i,1],pos_des[i,2]]), zeroadding=None)
                xnew_desiredy, tracedes_desiredy =interpolate_trace(xnew1, tracedes1y, point_online1, xnew2, tracedes2y, point_online2, np.array([0,pos_des[i,1],pos_des[i,2]]), zeroadding=None)
                xnew_desiredz, tracedes_desiredz =interpolate_trace(xnew1, tracedes1z, point_online1, xnew2, tracedes2z, point_online2, np.array([0,pos_des[i,1],pos_des[i,2]]), zeroadding=None)

                #print("langth traces:", len(txt0.T[0]), len(txt0.T[1]),len(txt1.T[0]),len(txt1.T[1]),len(txt2.T[0]), len(txt2.T[1]),len(txt3.T[0]),len(txt3.T[1]))
                #print("langth interpolated:", len(xnew1), len(tracedes1x),len(xnew2),len(tracedes2x),len(xnew_desiredx), len(tracedes_desiredx))

                #now, lets use timing solution
                t0=DesiredT0[i]
                ntbins=len(tracedes_desiredx)

                xnew_desiredx=np.linspace(t0+tmin,t0+tmax-tbinsize*3,ntbins)
                xnew_desiredy=xnew_desiredx
                xnew_desiredz=xnew_desiredx

                if(len(xnew_desiredx)!=ntbins):
                 print("warning! different sizes",ntbins,len(xnew_desiredx))
                 print(tmin,tmax,tbinsize,ntbins)

                if(round((xnew_desiredx[2]-xnew_desiredx[1]),3)!=round(tbinsize,3)):
                 print("warning! different tbin sizes",tbinsize,xnew_desiredx[2]-xnew_desiredx[1])
                 print(tmin,tmax,tbinsize,ntbins)


                if DISPLAY>2:

                    if(PLOTPAPER):
                     plt.rc('font', family='serif', size=15)

                    fig4 = plt.figure()
                    ax4 = fig4.add_subplot(1,2,2)
                    ax4.plot(txt0.T[0]/1000, txt0.T[2], label = "I")
                    ax4.plot(txt2.T[0]/1000, txt2.T[2], label = "II")
                    ax4.plot(txt3.T[0]/1000, txt3.T[2], label = "III")
                    ax4.plot(txt1.T[0]/1000, txt1.T[2], label = "IV")
                    ax4.plot(xnew1/1000, tracedes1y, linestyle='--',color='r',label = "I->IV")
                    ax4.plot(xnew2/1000, tracedes2y, linestyle='--', color='b',label = "II->III")
                    ax4.plot(xnew_desiredx/1000, tracedes_desiredy, linestyle='--',color='k', label = "Synthetized")
                    tmp=ax4.set(title="Signal Y component",xlabel='Time [$\mu$s]',ylabel='Amplitude')

                    plt.legend(loc= "best")

                    ax3 = fig4.add_subplot(1,2,1)

                    ax3.scatter(pos_sims[:,1], pos_sims[:,2], color='c',label = "Simulated")
                    ax3.scatter(pos_sims[points_I[0][0],1], pos_sims[points_I[0][0],2], marker ="D",label = "I")
                    ax3.scatter(pos_sims[points_II[0][0],1], pos_sims[points_II[0][0],2],marker ="D", label = "II")
                    ax3.scatter(pos_sims[points_III[0][0],1], pos_sims[points_III[0][0],2], marker ="D", label = "III")
                    ax3.scatter(pos_sims[points_IV[0][0],1], pos_sims[points_IV[0][0],2], marker ="D",label = "IV")
                    tmp=ax3.set(title="Desired Antenna:"+str(i),xlabel='vxB',ylabel='vx(vxB)')

                    #print("desired "+str(pos_des[i]))
                    #print("pointI "+str(points_I[0]))
                    #print("pointII "+str(points_II[0]))
                    #print("pointIII "+str(points_III[0]))
                    #print("pointIV "+str(points_IV[0]))

                    ax3.scatter(point_online1[1], point_online1[2], marker ="x",color='r',label = "I->IV")
                    ax3.scatter(point_online2[1], point_online2[2], marker ="x",color='b',label = "II->III")
                    ax3.scatter(pos_des[i,1], pos_des[i,2], color= 'k',label = "Synthetized")

                    ax3.scatter(0, 0, marker ="x",  label = "core")
                    ax3.plot([0, pos_des[i,1]], [0, pos_des[i,2]])

                    for j in range(0,len(pos_sims[:,1])):
                        ax3.annotate(str(j), ((pos_sims[j,1], pos_sims[j,2])),fontsize=10)

                    plt.legend(loc= "best")
                    #plt.show()

                if DEVELOPMENT and DISPLAY>2:

                   AntennaID=hdf5io.GetAntennaID(CurrentAntennaInfo,160+i)
                   if(tracetype=='efield'):
                      txtdes=hdf5io.GetAntennaEfield(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                   elif(tracetype=='voltage'):
                      txtdes=hdf5io.GetAntennaVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                   elif(tracetype=='filteredvoltage'):
                      txtdes=hdf5io.GetAntennaFilteredVoltage(InputFilename,CurrentEventName,AntennaID,OutputFormat="numpy")
                   else:
                      print('You must specify either efield, voltage or filteredvoltage, bailing out')
                      return 0

                   fig5, ((ax5, ax6), (ax7, ax8), (ax9, ax10))  = plt.subplots(3, 2, sharey='row', sharex='col')

                   ax5.plot(xnew_desiredx/1000, tracedes_desiredx, linestyle='--',color='g', label = "Synthetized")
                   ax5.plot(txtdes.T[0]/1000, txtdes.T[1], label = "Simulation")
                   tmp=ax5.set(ylabel='Signal in X')
                   #ax5.legend(loc='lower right')

                   if(tracetype=='efield'):
                     difference=tracedes_desiredx-txtdes.T[1,0:-1]
                     time=txtdes.T[0,0:-1]/1000
                   else:
                     difference=tracedes_desiredx-txtdes.T[1]
                     time=txtdes.T[0]/1000

                   ax6.plot(time,difference, linestyle='--',color='g', label = "Difference")
                   #ax6.legend(loc='lower right')

                   ax7.plot(xnew_desiredy/1000, tracedes_desiredy, linestyle='--',color='g', label = "Synthetized")
                   ax7.plot(txtdes.T[0]/1000, txtdes.T[2], label = "Simulation")
                   tmp=ax7.set(ylabel='Signal in Y')
                   #ax7.legend(loc='lower right')

                   if(tracetype=='efield'):
                     difference=tracedes_desiredy-txtdes.T[2,0:-1]
                     time=txtdes.T[0,0:-1]/1000
                   else:
                     difference=tracedes_desiredy-txtdes.T[2]
                     time=txtdes.T[0]/1000

                   ax8.plot(time, difference, linestyle='--',color='g', label = "Difference")
                   #ax8.legend(loc='lower right')

                   ax9.plot(xnew_desiredz/1000, tracedes_desiredz, linestyle='--',color='g', label = "Synthetized")
                   ax9.plot(txtdes.T[0]/1000, txtdes.T[3], label = "Simulation")
                   tmp=ax9.set(ylabel='Signal in Z',xlabel='Time [$\mu$s]')
                   ax9.legend(loc='lower right')

                   if(tracetype=='efield'):
                     difference=tracedes_desiredz-txtdes.T[3,0:-1]
                     time=txtdes.T[0,0:-1]/1000
                   else:
                     difference=tracedes_desiredz-txtdes.T[3]
                     time=txtdes.T[0]/1000

                   #ax10 = fig5.add_subplot(3,2,6)
                   ax10.plot(time, difference, linestyle='--',color='g', label = "Difference")
                   tmp=ax10.set(xlabel='Time [$\mu$s]')
                   ax10.legend(loc='lower right')

                if DISPLAY>2:
                  plt.show()

                # ------------------

                # Save as textfile
                #print("Interpolated trace stord as ",split(desired)[0]+ '/a'+str(ind[i])+'.trace')
                #os.system('rm ' + split(desired)[0] + '/*.trace')
                #FILE = open(split(desired)[0]+ '/a'+str(ind[i])+'.trace', "w+" )

                #uncomment this if you want the old style output
                #FILE = open(directory+ '/Test/a'+str(ind[i])+'.trace', "w+" )
                #for j in range( 0, len(xnew_desiredx) ):
                #        print("%3.2f %1.5e %1.5e %1.5e" % (xnew_desiredx[j], tracedes_desiredx[j], tracedes_desiredy[j], tracedes_desiredz[j]), end='\n', file=FILE)
                #FILE.close()

                if(tracetype=='efield'):
                    efield=np.column_stack((xnew_desiredx,tracedes_desiredx,tracedes_desiredy,tracedes_desiredz))
                    EfieldTable=hdf5io.CreateEfieldTable(efield, CurrentEventName, CurrentEventNumber , DesiredIds[i], i, "Interpolated", info={})
                    hdf5io.SaveEfieldTable(OutputFilename,CurrentEventName,str(DesiredIds[i]),EfieldTable)
                elif(tracetype=='voltage'):
                    voltage=np.column_stack((xnew_desiredx,tracedes_desiredx,tracedes_desiredy,tracedes_desiredz))
                    VoltageTable=hdf5io.CreateVoltageTable(voltage, CurrentEventName, CurrentEventNumber , DesiredIds[i], i, "Interpolated", info={})
                    hdf5io.SaveVoltageTable(OutputFilename,CurrentEventName,str(DesiredIds[i]),VoltageTable)
                elif(tracetype=='filteredvoltage'):
                    filteredvoltage=np.column_stack((xnew_desiredx,tracedes_desiredx,tracedes_desiredy,tracedes_desiredz))
                    VoltageTable=hdf5io.CreateVoltageTable(filteredvoltage, CurrentEventName, CurrentEventNumber , DesiredIds[i], i, "Interpolated", info={})
                    hdf5io.SaveFilteredVoltageTable(OutputFilename,CurrentEventName,str(DesiredIds[i]),VoltageTable)


                #delete after iterate
                del points_I, points_II, points_III, points_IV


    #now, lets remove the antennas fromthe index
    #for i in remove_antenna:
    DesiredIds=np.delete(DesiredIds,remove_antenna)
    DesiredAntx=np.delete(DesiredAntx,remove_antenna)
    DesiredAnty=np.delete(DesiredAnty,remove_antenna)
    DesiredAntz=np.delete(DesiredAntz,remove_antenna)
    DesiredSlopeA=np.delete(DesiredSlopeA,remove_antenna)
    DesiredSlopeB=np.delete(DesiredSlopeB,remove_antenna)
    DesiredT0=np.delete(DesiredT0,remove_antenna)
    #CreateAntennaInfo(IDs, antx, anty, antz, antt, slopeA, slopeB, AntennaInfoMeta, P2Pefield=None,P2Pvoltage=None,P2Pfiltered=None,HilbertPeak=None,HilbertPeakTime=None):
    DesiredAntennaInfo=hdf5io.CreateAntennaInfo(DesiredIds, DesiredAntx, DesiredAnty, DesiredAntz, DesiredT0, DesiredSlopeA, DesiredSlopeB, DesiredAntennaInfoMeta)
    hdf5io.SaveAntennaInfo(OutputFilename,DesiredAntennaInfo,CurrentEventName)

    #now aim at the point where all the antennas where interpolated and saved to file. Now i will calulate the peak to peak and hilbert envelope peak and time
    #this is done after everything was computed.
    if(usetrace=="all"):
      print("Computing P2P for "+str(OutputFilename))

      OutAntennaInfo=hdf5io.GetAntennaInfo(OutputFilename,CurrentEventName)
      OutIDs=hdf5io.GetAntIDFromAntennaInfo(OutAntennaInfo)

      p2pE=hdf5io.get_p2p_hdf5(OutputFilename,usetrace='efield')
      p2pV=hdf5io.get_p2p_hdf5(OutputFilename,usetrace='voltage')
      p2pFV=hdf5io.get_p2p_hdf5(OutputFilename,usetrace='filteredvoltage')

      peaktimeE, peakE=hdf5io.get_peak_time_hilbert_hdf5(OutputFilename,usetrace='efield')
      peaktimeV, peakV=hdf5io.get_peak_time_hilbert_hdf5(OutputFilename,usetrace='voltage')
      peaktimeFV, peakFV=hdf5io.get_peak_time_hilbert_hdf5(OutputFilename,usetrace='filteredvoltage')

      AntennaP2PInfo=hdf5io.CreateAntennaP2PInfo(OutIDs, DesiredAntennaInfoMeta, P2Pefield=p2pE,P2Pvoltage=p2pV,P2Pfiltered=p2pFV,HilbertPeakE=peakE,HilbertPeakV=peakV,HilbertPeakFV=peakFV,HilbertPeakTimeE=peaktimeE,HilbertPeakTimeV=peaktimeV,HilbertPeakTimeFV=peaktimeFV)
      hdf5io.SaveAntennaP2PInfo(OutputFilename,AntennaP2PInfo,CurrentEventName)


#-------------------------------------------------------------------

def interpol_check_hdf5(InputFilename, positions, new_pos, p2pE, InterpolMethod,usetrace='efield', DISPLAY=False):
    '''
    Interpolates the signal peak-to-peak electric field at new antenna positions
    Check that the interpolation efficiency at 6 antenna positions available in each shower file

    Parameters:
    InputFilename: str
        HDF5File
    positions: numpy array
        x, y, z coordinates of the antennas in the simulation (not used in trace interpolation method)
    new_pos: numpy array
        x, y, z coordinates of the antennas in new layout (at 6 check points)
    p2pE: numpy array
        [p2p_Ex, p2p_Ey, p2p_Ez, p2p_total]: peak-to-peak electric fields along x, y, z, and norm

    InterpolMethod: str
        interpolation method
        'lin' = linear interpolation from scipy.interpolate
        'rbf' = radial interpolation from scipy.interpolate
        'trace' = interpolation of signal traces:
            generates new interpolated trace files in path/Test/ directory

    DISPLAY: boolean
        if TRUE: 2D maps of peak-to-peak electric field
            at original and interpolated antennas are displayed

    Output:
    interp_err: numpy arrays
        interpolation error at each antenna (interpolated - original)/original
    p2p_total_new: numpy array
        peak-to-peak electric field at new antenna positions

    '''


    # interpolate (check rbf)
    logging.debug('interpol_check:Interpolating...'+str(usetrace))
    #print('Interpolating...'+path)

    number_ant = 160
    icheck = np.mgrid[160:176:1]

    myx_pos = positions[0,0:number_ant-1]
    myy_pos = positions[1,0:number_ant-1]
    myz_pos = positions[2,0:number_ant-1]
    mypositions = np.stack((myx_pos, myy_pos, myz_pos), axis=0)
    myp2p_total = p2pE[3,0:number_ant-1]

    from trace_interpol_hdf5 import do_interpolation_hdf5
    OutputFilename = InputFilename + '.Interpolated.'+str(usetrace)+'.hdf5'

    #do_interpolation(AntPath,new_pos,mypositions,Zenith,Azimuth,phigeo=147.43, thetageo=0.72, shower_core=np.array([0,0,2900]), DISPLAY=False)
    do_interpolation_hdf5(new_pos, InputFilename, OutputFilename, antennamin=0, antennamax=159, EventNumber=0, DISPLAY=DISPLAY, usetrace=usetrace)

    #NewAntNum = size(new_pos)
    #NewAntNum, NewAntPos, NewAntID = get_antenna_pos_zhaires(NewAntPath)
    #NewP2pE = get_p2p(path+"/Test",NewAntNum)

    NewP2pE = hdf5io.get_p2p_hdf5(OutputFilename,antennamax=15,antennamin=0,usetrace=usetrace)

    p2p_total_new = NewP2pE[3,:]
    p2p_x_new = NewP2pE[0,:]
    p2p_y_new = NewP2pE[1,:]
    p2p_z_new = NewP2pE[2,:]

    # checking the interpolation efficiency
    interp_err = abs(p2p_total_new-p2pE[3,icheck])/p2pE[3,icheck]
    interp_errx = abs(p2p_x_new-p2pE[0,icheck])/p2pE[0,icheck]
    interp_erry = abs(p2p_y_new-p2pE[1,icheck])/p2pE[1,icheck]
    interp_errz = abs(p2p_z_new-p2pE[2,icheck])/p2pE[2,icheck]

    #print(np.shape(p2p_total_new))
    #print(np.shape(p2pE[3,icheck]))
    #print(p2pE[3,icheck])
    #print("interp_err = #{}".format(interp_err))


    if (DISPLAY and InterpolMethod!='trace'):
        logging.debug('interpol_check:Plotting...')

        ##### Plot 2d figures of total peak amplitude in positions along North-South and East-West
        fig1 = plt.figure(10,figsize=(5,7), dpi=100, facecolor='w', edgecolor='k')


        ax1=fig1.add_subplot(211)
        name = 'total'
        plt.title(name)
        ax1.set_xlabel('positions along NS (m)')
        ax1.set_ylabel('positions along EW (m)')
        col1=ax1.scatter(positions[0,:],positions[1,:], c=p2pE[3,:],  vmin=min(myp2p_total), vmax=max(myp2p_total),  marker='o', cmap=cm.gnuplot2_r)
        plt.xlim((min(mypositions[0,:]),max(mypositions[0,:])))
        plt.ylim((min(mypositions[1,:]),max(mypositions[1,:])))
        plt.colorbar(col1)
        plt.tight_layout()


        ax2=fig1.add_subplot(212)
        name = 'total interpolated'
        plt.title(name)
        ax2.set_xlabel('x (m)')
        ax2.set_ylabel('y (m)')
        col2=ax2.scatter(new_pos[0,:],new_pos[1,:], c=p2p_total_new,  vmin=np.min(myp2p_total), vmax=np.max(myp2p_total),  marker='o', cmap=cm.gnuplot2_r)
        plt.xlim((min(mypositions[0,:]),max(mypositions[0,:])))
        plt.ylim((min(mypositions[1,:]),max(mypositions[1,:])))
        plt.colorbar(col2)
        plt.tight_layout()


        plt.show(block=False)


        if (interp_err.min() < 1.e-9):
            fig2 = plt.figure(figsize=(5,7), dpi=100, facecolor='w', edgecolor='k')


            ax1=fig2.add_subplot(211)
            name = 'total'
            plt.title(name)
            ax1.set_xlabel('positions along NS (m)')
            ax1.set_ylabel('positions along EW (m)')
            col1=ax1.scatter(positions[0,:],positions[1,:], c=p2pE[3,:],  vmin=min(myp2p_total), vmax=max(myp2p_total),  marker='o', cmap=cm.gnuplot2_r)
            plt.xlim((min(mypositions[0,:]),max(mypositions[0,:])))
            plt.ylim((min(mypositions[1,:]),max(mypositions[1,:])))
            plt.colorbar(col1)
            plt.tight_layout()


            ax2=fig1.add_subplot(212)
            name = 'total interpolated'
            plt.title(name)
            ax2.set_xlabel('x (m)')
            ax2.set_ylabel('y (m)')
            col2=ax2.scatter(new_pos[0,:],new_pos[1,:], c=p2p_total_new,  vmin=np.min(myp2p_total), vmax=np.max(myp2p_total),  marker='o', cmap=cm.gnuplot2_r)
            plt.xlim((min(mypositions[0,:]),max(mypositions[0,:])))
            plt.ylim((min(mypositions[1,:]),max(mypositions[1,:])))
            plt.colorbar(col2)
            plt.tight_layout()


            plt.show(block=False)



    return interp_err, p2p_total_new, interp_errx, p2p_x_new, interp_erry, p2p_y_new, interp_errz, p2p_z_new





def main():
    if ( len(sys.argv)<1 ):
        print("""
            Example on how to do interpolate a signal
                -- read in list of desired poistion
                -- read in already simulated arrazs
                -- find neigbours and perform interpolation
                -- save interpolated trace

            Usage: python3 interpolate.py <path>
            Example: python3 interpolate.py <path>

            path: Filename and path of the input hdf5file
        """)
        sys.exit(0)

    hdf5file = sys.argv[1]

    # path to list of desied antenna positions, traces will be stored in that corresponding folder
    #desired  = sys.argv[1]
    #desired=np.array([[ 100., 0., 2900.],[ 0., 100., 2900.]])
    desired=np.loadtxt("/home/mjtueros/GRAND/GP300/GridShape/Stshp_XmaxLibrary_0.1995_85.22_0_Iron_23/Test/new_antpos.dat",usecols=(2,3,4))

    OutputFilename=split(hdf5file)[0]+"/InterpolatedAntennas.hdf5"
    print(OutputFilename)
    # call the interpolation: Angles of magnetic field and shower core information needed, but set to default values
    #do_interpolation(desired,hdf5file, zenith, azimuth, phigeo=147.43, thetageo=0.72, shower_core=np.array([0,0,2900]), DISPLAY=False)
    do_interpolation_hdf5(desired, hdf5file, OutputFilename, antennamin=0, antennamax=159, EventNumber=0, DISPLAY=False,usetrace='efield')

if __name__== "__main__":
  main()

