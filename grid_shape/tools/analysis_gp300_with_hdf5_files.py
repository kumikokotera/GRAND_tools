import hdf5fileinout as hd


import matplotlib.pyplot as plt
import numpy as np
import os
import json
from grid_shape import diff_spec as diff_spec
from grid_shape import grids as grids
from grid_shape import utils_analysis as ua
 
radius = 125

base = grids.create_grid_univ("trihex", 125, input_n_ring=50)


gp300radius = (1 + 90 * 1.5) * 2 / np.sqrt(3) * 125 + radius #kewen original array is 90 rings hexagon with step 125
#print("hexagon radius", gp300radius,radius)
found = False
while(not found):
    offset = grids.get_offset_in_grid("hexhex", gp300radius)
    # offset=[-12700,11900] #this will give 00 ofset to the north, 2000 to the east (to put in NE corner: -12700,11900)
    # move coordinates to ftp
    # change x and y to get to the hexagon in grand coordinates
    offset = [offset[1], offset[0]]
    xftp = offset[0] + 5.0 * 1000.0
    yftp = offset[1] + 7.7428 * 1000.0
    #outside of radius?
    #x^2+y^2=25
    if (xftp**2 + yftp**2 > (25000 + radius)**2):
    #print(offset,xftp,yftp,"is outside the 25km limit")
    #outsidePoints1=np.vstack((outsidePoints1,offset))             
        continue
    #superior N
    #y>-4.0730*x+64.4106 (x,y en km) + c ->  y+4.0730x-64410.6 > c  (x,y en metros)              
    if (yftp + 4.073 * xftp - 64410.6 > radius / np.cos(np.arctan(4.073))):
    #print(offset,xftp,yftp,"is above the northern limit",yftp+4.073*xftp-64410.6)
    #outsidePoints2=np.vstack((outsidePoints2,offset))              
        continue
    #inferior SW
    #y<-1.765*x+0.8825 - c (x,y en km)->  y<-1.765*x+882.5 -c 
    #inferior SE 
    #y>2.1471*x-1.0735 + c (x,y en km) ->  y>2.1471*x - 1073.5 + c
    if (yftp + 1.765 * xftp - 882.5 < -radius/np.cos(np.arctan(1.765)) and yftp - 2.1471 * xftp + 1073.5 > radius/np.cos(np.arctan(2.1471))):
    #print(offset,xftp,yftp,"is below the suthern limit")
    #outsidePoints3=np.vstack((outsidePoints3,offset))              
        continue
    found = True
    print("using grid point offset to generate random core",offset)  

    