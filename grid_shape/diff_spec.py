import os
import numpy as np



def tale_diff_flux(E): # E en eV
    """
    spectrum data from TALE https://arxiv.org/pdf/1803.01288.pdf
    input: energy in eV

    output: J(E) in eV^-1 m^-2 s^-1 sr^-1 
    """

  
    bp1 = 16.22
    bp2 = 17.04

    gamma3 = -3.19
    J3E3 = 3.504 * 1E24 ;
    E3 = 10**17.05
    A3 = (J3E3/(E3**3))/E3**gamma3;

    gamma2=-2.92
    J2E3 = 3.504*1E24
    E2 = 10**17.05
    A2 = (J2E3/(E2**3))/(E2**gamma2);
   
    gamma1 = -3.12;
    J1E3 = 3.405*1E24;
    E1 = 10**15.7;
    A1 = (J1E3/(E1**3)) / (E1**gamma1);

    E_log = np.log10(E)

    if E_log >= bp2:
        return A3 * E**gamma3
    elif (E_log >= bp1) * (E_log < bp2):
        return A2 * E**gamma2
    elif (E_log < bp1):
        return A1 * E**gamma1
    else:
        return -1    

