# Author: Nicole Farias
# Date: July 2021
# Description: This code calculates the transmission through a surface (located in the z plane)
# from the directivity data provided by HFSS
# It assumes that any gain from theta < 90 is reflected power, whereas transmitted power comes from 
# gains with theta > 90


#Importing relevant libraries
import csv
import numpy as np
from matplotlib import pyplot as plt
import math

pi = math.pi

# Main function
def main():
    
    print("\n\n-------------\nWelcome. Let's see if stuff fits!")
    
    R = 16000 # um -- radius of lenslet
    
    #metamaterial parameters
    TH = 700 # top height
    BH = 545 # bottom height
    
    TD = 160 # top diameter 
    BD = 120 # bottom diameter
    s = 10 # spacing between holes
    p = TD + s # pitch
    
    # number of holes:
    N = 2*np.pi*R/p
    
    # total hole depth
    H = TH+BH
    
    # Available circumference:
    
    s0 = get_s_r(R, TD, N)
    
    s1 = get_s_r(R-TH, TD, N)
    
    s2 = get_s_r(R-H, BD, N)
    
    print("total height: %.2f" %(H))
    print("S0: %.2f,  S1: %.2f  , S2: %.2f " %(s0, s1, s2))
    
    
    # F1: 10.00,  S1: 2.56  , S2: 36.77]v\
    # F2: 10.00,  S1: 1.50  , S2: 34.59
    # F3: S0: 10.00,  S1: 1.78  , S2: 35.12
    
    print("\nEnd of program. \n------------\n\n")
    
def get_s_r(r, d, N):
    # radius from lenslet center
    # diameter of cylinder
    c = circ(r) # the circumference at radius r
    p_avail = c/N # pitch available
    s = p_avail-d # space between holes
    return s
    
def circ(r):
    # r: radius
    # get the circumference
    return 2*np.pi*r
    
    

if __name__ == "__main__":
    main()