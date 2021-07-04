#This program does a simple fit of transmission data to get the thickness and index of refraction of anti-reflection coatings

import math
import numpy as np
import cmath
import matplotlib.pyplot as plt

def main():

    print("Welcome! This program calculates the index of refraction and thickness of anti-reflection coatings \nby doing a simple fit of transmission data. \nYou must provide the data file and the ranges of interest.")

    # DATA FILE LOCATION
    filepath = 'C:/Users/nicol/Documents/00Research/PythonCode/AR simulation/'
    filename = 'ctthin_transmission.txt'
    #Read data file
    freq_array_raw, T_array_raw = ReadTransmissionData(filepath, filename)
    #print(freq_array_raw)

    #CONSTANTS
    #indices of refracion
    n_vacuum = 1 # index of refraction of vacuum
    n_substrate = 3.11 # index of refraction of substrate/lens. Typically known. 3.11 for alumina, 3.4 for silicon
    lt_substrate = 0.00001 # loss tangent of substrate/lens

    #RANGE OF EXPECTED INDEX OF REFRACTION
    n_range = np.linspace(1.4, 3.2, num=20)
    #RANGE OF EXPECTED LOSS TANGENT
    lt_range = np.linspace(0, 1e-4, num=20)
    

    # angle of incidence
    theta0 = 0 

    # thickness of samples:
    t_substrate = 0.25 * 2.54/100 # [m] thickness of alumina witness sample
    t_ar = np.array([0.017, 0.0105]) * 2.54/100 # [m] thickness of AR layers (provided by Oliver)


    # Initial guess of index of refraction of ar coating
    n_ar = [math.sqrt(2), math.sqrt(5.2)]

    # SETUP LAYERS
    optics_type = "lens" # either "lens" or "lenslet"
    
    if optics_type == "lens":
        n_ar_flipped = np.flipud(n_ar)
        t_ar_flipped = np.flipud(t_ar)
        n_array = [n_vacuum, n_ar, n_substrate, n_ar_flipped, n_vacuum]
        t_array = [0, t_ar, t_substrate, t_ar_flipped]

    print(n_array)
    print(t_array)



    print("\nThe end.\n\n")

def ReadTransmissionData(filepath: str, filename: str):
    fullname = filepath + filename
    with open(fullname) as file:
        data = file.readlines()
        freq_array = np.zeros(len(data)) # array to store frequency array
        T_array = np.zeros(len(data)) # array to store transmission data                                 
        f = 0
        for line in data:
            freqS, TS, temp = (line.split()) #strings third column is useless
            freq_array[f] = float(freqS) * 10**9 # convert from Ghz to Hz
            T_array[f] = float(TS)
            f = f + 1 # a counter
                                                                                         
    return freq_array, T_array


   




if __name__ == "__main__":
    main()

