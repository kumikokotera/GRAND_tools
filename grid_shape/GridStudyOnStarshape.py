'''
Count triggered antennas for one shower, for a set of layouts
Performs an interpolation of the traces of the peak-to-peak electric field/voltage/filtered voltage
at new layout antenna positions

Test with $PYHONINTERPRETER GridStudyOnStarshape.py /home/mjtueros/AiresRepository/DiscreteLibrary/TopoStshpOutbox/Stshp_Proton_3.98_86.5_0.0_9/Stshp_Proton_3.98_86.5_0.0_9.hdf5 .

'''
#!/usr/bin/env python
import sys
import os
import logging
import glob
import numpy as np
ZHAIRESRUNNER=os.environ["ZHAIRESRUNNER"]
sys.path.append(ZHAIRESRUNNER + "/ConicalInterpolator") #so that it knows where to find things
from StarshapeInterpolation import do_interpolation_hdf5
ZHAIRESPYTHON=os.environ["ZHAIRESPYTHON"]
sys.path.append(ZHAIRESPYTHON)
import hdf5fileinout as hdf5io
sys.path.append("/home/mjtueros/GRAND/GP300/GridShapeKotera/GRAND_tools/grid_shape")
import grids


def grid_study_on_starshape_hdf5(InputDir, OutputDir, configurations,  usetrace='efield'):

    DISPLAY=False

    #this could be done more elegantly by using the status file if we delete Aires.status from the standard, or we open any hdf5 file looking for the JobName...
    InputFilename=glob.glob(InputDir+"/*.hdf5")

    for i in range(0,len(InputFilename)):
      if "Interpolated" not in InputFilename[i] and "NoTraces" not in InputFilename[i]:
        InputFileName=InputFilename[i]
        break

    #get basi showerinfo we will output to the showerinfo file
    CurrentRunInfo=hdf5io.GetRunInfo(InputFileName)
    CurrentEventName=hdf5io.GetEventName(CurrentRunInfo,0) #using the first event of each file (there is only one for now)
    InputJobName=CurrentEventName
    CurrentAntennaInfo=hdf5io.GetAntennaInfo(InputFileName,CurrentEventName)
    Zenith=hdf5io.GetEventZenith(CurrentRunInfo,0)
    Azimuth=hdf5io.GetEventAzimuth(CurrentRunInfo,0)
    Primary=hdf5io.GetEventPrimary(CurrentRunInfo,0)
    Energy=hdf5io.GetEventEnergy(CurrentRunInfo,0)
    XmaxDistance=hdf5io.GetEventXmaxDistance(CurrentRunInfo,0)
    SlantXmax=hdf5io.GetEventSlantXmax(CurrentRunInfo,0)
    HadronicModel=hdf5io.GetEventHadronicModel(CurrentRunInfo,0)
    RandomAzimuth=np.random.uniform(0,180)

    print(RandomAzimuth) #in this version all will share the same azimuth, which is more desirable

    for case in configurations:
        gridtype=case[0]
        step=case[1]

        NewPos,RandomCore = grids.create_grid_univ(gridtype,step,RandomAzimuth,do_offset=True,DISPLAY=DISPLAY)
        NewPos=NewPos.T

        print("NewPos",NewPos,len(NewPos[:,0]))
        print("RandomCore",RandomCore)

        #Now The interpolation requests the T0 of the antennas. I need to work on that.
        #for now, this is dummy

        Dummy= np.random.rand(len(NewPos[:,0]),4)
        Dummy[:,1:4]=NewPos
        NewPos=Dummy


        OutputFileName= OutputDir+"/"+InputJobName + ".Interpolated."+gridtype+"_"+str(step)+"_"+str(usetrace)+".hdf5"

        print("about to do interpolation output to:"+str(OutputFileName))

        do_interpolation_hdf5(NewPos, InputFileName, OutputFileName, antennamin=0, antennamax=159, EventNumber=0,DISPLAY=DISPLAY, usetrace=usetrace)

        if(usetrace=="all"):
          print("usetrace is all, looping over all trace types")
          usetracelist=["efield","voltage","filteredvoltage"]
        else:
          usetracelist=[str(usetrace)]

        for tracetype in usetracelist:

          p2p_out = hdf5io.get_p2p_hdf5(OutputFileName,usetrace=tracetype)

          #now i have to complete with 0 the antennas of the array that where not interpolated, so that the output is easier to use
          OutRunInfo=hdf5io.GetRunInfo(OutputFileName)
          OutEventName=hdf5io.GetEventName(OutRunInfo,0)
          OutAntennaInfo=hdf5io.GetAntennaInfo(OutputFileName,OutEventName)
          OutIDs=hdf5io.GetAntIDFromAntennaInfo(OutAntennaInfo)

          p2p_total_new=np.zeros((4,len(NewPos)))

          index=0
          for i in OutIDs:
            p2p_total_new[0,int(i)]=p2p_out[0,index]
            p2p_total_new[1,int(i)]=p2p_out[1,index]
            p2p_total_new[2,int(i)]=p2p_out[2,index]
            p2p_total_new[3,int(i)]=p2p_out[3,index]
            index=index+1

          # write antenna trigger information file (this would be cool to have it also inside the interpolated hdf5)
          P2PFile = OutputDir+"/"+InputJobName + ".Interpolated."+gridtype+"_"+str(step)+"_"+str(tracetype)+".P2Pdat"
          np.savetxt(P2PFile, p2p_total_new,fmt='%1.4e')

          FILE = open(P2PFile+str(".showerinfo"),"w" )
          print("%s %s %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %s" % (InputJobName,Primary,Energy,Zenith,Azimuth,XmaxDistance,SlantXmax,RandomCore[0],RandomCore[1],RandomAzimuth,HadronicModel), file=FILE)
          FILE.close()

        print("done with "+str(case))

    print("done with all configurations")
    #here i could clean the files i dont want to keep:
    #in this case, i will delete all the hdf5 files, becouse they take too much space and i dont have any use for them for now
    for name in glob.glob(OutputDir+"/"+InputJobName + ".Interpolated.*.hdf5"):
      os.remove(name)



def main():

    if ( len(sys.argv)<3 ):
        print("""

            Usage: python3  inputdirectory outputdirectory

        """)
        sys.exit(0)

    InputDir = sys.argv[1]
    OutputDir = sys.argv[2]

    # create rectangular antenna layouts with various steps
    #rectStep = [137.5, 275, 550,825, 1100] # in m
    rectStep = [137.5] # in m

    # create hexagonal antenna layouts with various radii
    #hexStep = [125, 250, 500, 750, 1000] # radius of hexagon in m
    hexStep = [125] # radius of hexagon in m

    configurations=[["trihex",125],["hexhex",125],["rect",137.5],
                    ["trihex",250],["hexhex",250],["rect",275],
                    ["trihex",500],["hexhex",500],["rect",500],
                    ["trihex",750],["hexhex",750],["rect",825],
                    ["trihex",1000],["hexhex",1000],["rect",1100]]

    #configurations=[["trihex",125],["hexhex",125],["rect",137.5]]


    grid_study_on_starshape_hdf5(InputDir, OutputDir, configurations, usetrace='all')


print(__name__)
if __name__ == '__main__':
    main()
