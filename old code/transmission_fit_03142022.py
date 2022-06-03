# This code fits the ar coating parameters of a given transmissiond data
# Author: Nicole Farias
# Last updated: 09 09 2021

# TODOs:
# 1 - allow all lts to be inserted from json file 
# 1.1 - setup the lts array in the "Setup Layers" function
# 2 - allow for 2 layers of coatings (see Fit function)

import math
import numpy as np
import cmath
import matplotlib.pyplot as plt
import CharlesHill_PB2b_AR_optimize as Charlie
from transmission_data_analysis import *
import json

# global constants
N_VACUUM = 1

def main():
    
    print("\n\nWelcome! This program calculates the index of refraction and thickness of anti-reflection coatings \nby doing a simple fit with many for loops :P ")
    
    # fit parameters
    tol = 0.2 # 20%
    N = 10 # number of data points to be analyzed per fit.
    
    # DATA FILE LOCATION
    filepath = "C:/Users/nicol/Documents/00Research/Data/Thermal_Sprayed_Samples/2021_august_CTS_Hayward_first_batch/"
    filename = "120a_G_trimmed"
    ftype = ".csv"
    f_loc = filepath + filename + ftype
    json_loc = filepath + filename + '.txt'
    
    # JASON FILE LOCATION
    with open(json_loc) as json_file:
        D = json.load(json_file) # data
        print(D)
        #print(data['filename'])
    
    
    
    # READ THE DATA FILE
    if ftype == ".csv":
        freq_array_raw, T_array_raw, var_array, title = read_csv_file(f_loc, D['freq_column'],  D['data_column'], D['var_column'])
        print("Read data sucessfully.")
        
    else:
        print("Data file type not recognized.")

    
    # UNCOMMENT THIS IF USING SIMULATION DATA
    #split data if necessary (i.e., in simulations we might be looking at different geometries, with a varying radius or something like that)
    #if D['var_column'] >= 0:
    #    freq_array_raw, T_array_raw, var_values  = split_data(var_array, freq_array_raw, T_array_raw) # I have to check where this function is located!
    #    #select one of the columns
    #    T_array_raw = T_array_raw[:, var_selection]
    #    var_label = str(var_values[var_selection])
    #else:
    #    var_values = []

    # Transform frequency into Hz
    freq_array_raw = freq_array_raw * D['f_convert']  # first, convert everything into Hz
    # TRIM THE DATA IF NECESSARY (if we only want to look at a certain frequency band)
    freq_array, T_array_data = SelectFrequencyBand(D['freq_min'], D['freq_max'], freq_array_raw, T_array_raw) # this is used in case we want to zoom in at a specific frequency
    print("Range of frequency: ", int(min(freq_array/1e9)), " - " , int(max(freq_array/1e9)) , " GHz")

    # create the arrays with AR characteristics
    t_ar = np.array([]) # thickness
    n_ar = np.array([]) # refraction index
    if D['t_ar_1'] >0: # if there is an AR layer
        t_ar = np.append(t_ar, D['t_ar_1'])
        n_ar = np.append(n_ar, D['n_ar_1'])
    if D['t_ar_2'] >0: # if there is a second AR layer
        t_ar = np.append(t_ar, D['t_ar_2'])
        n_ar = np.append(n_ar, D['n_ar_2'])


    #EVALUATE INITIAL GUESS
    n_array, t_array = SetupLayers(D['opt_type'], n_ar, N_VACUUM, D['n_substrate'], D['t_substrate'], t_ar)
    lt_array = [ 0.001500, 0.02, 0]
    #n_array = [1 , 3.11, 2.378 ,1]
    #t_array = [ 0, 0.00635, 0.00049073]
    #lt_array = [0.02, 0.001500, 0]
    Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, t_array, lt_array, freq_array)
    Ts_array_guess = 1 - Refl_s
    guess_error = CalculateError(T_array_data, Ts_array_guess)
    print("\n--------\nInitial guess: ")
    print("n_array: ", n_array)
    print("t_array: ", t_array)
    print("lt_array: ", lt_array)
    print("Error: ", guess_error)
    print("\n ------- \n ")
    #PLOT THE INITIAL GUESS
    plt.figure(1)
    plt.plot(freq_array/1e9, T_array_data, label="Data")
    plt.plot(freq_array/1e9, Ts_array_guess, label="Guess")
    plt.legend()
    plt.ylabel("Transmission Power")
    plt.xlabel("Frequency [GHz]")
    plt.title(D['data_label'])
    plt.show()

    # CALCULATE THE BEST FIT
    theta0 = np.pi/180*D['inc_angle'] # the incidence angle
    Ts_array_fit, n_array_fit, lt_array_fit, t_array_fit, error = FitDielectricConstant( T_array_data, freq_array, D, t_ar, n_ar, theta0, N, tol)

    # PLOT FIT AND DATA
    plt.figure(2)
    plt.plot(freq_array/1e9, T_array_data, label="Data")
    plt.plot(freq_array/1e9, Ts_array_fit, label="Fit")
    plt.plot(freq_array/1e9, Ts_array_guess, label="Guess")
    plt.legend()
    plt.ylabel("Transmission Power")
    plt.xlabel("Frequency [GHz]")
    plt.title(D['data_label'])
    
    
    # PRINT RESULTS
    print("\n----\nFit Results")
    print("N array: ", n_array_fit)
    print("lt_array_fit", lt_array_fit)
    print("t_array_fit: ", t_array_fit)
    print("error: ", error)
    
    plt.show()
    print("The end")
##################################################################

def FitDielectricConstant( T_array_data, freq_array, D, t_ar, n_ar, theta0:float, N, tol):
    """ Calculate a simple best fit of the dielectric constant
        ONLY WORKING FOR SINGLE LAYER AR
    Args:
        T_array_data ([type]): array containing transmission data
        freq_array ([type]): array containing frequency data
        n_array_0 ([type]): initial guess of indices of refraction
        lt_array_0 ([type]): initial guess of loss tangent
        d_array ([type]): array containing thickness of materials
        theta0 (float): angle of incidence
        N : number of trials per parameter
        tol : range of acceptable parameters (maybe change the name of this variable?)

    Returns:
        [type]: [description]
    """
    
    # initialize the error value with something large
    min_error = 1e5    
    n_coatings = len(t_ar)
    
    # SETUP THE RANGES ARRAYS
    #SUBSTRATE N AND T
    n_substrate = D['n_substrate']
    t_substrate = D['t_substrate']
    n_substrate_array = np.linspace( (1-tol) * n_substrate, (1+tol)*n_substrate, N)
    t_substrate_array = np.linspace( (1-tol) * t_substrate, (1+tol)*t_substrate, N)
    
    # AR COATING -- AT SOME POINT I NEED TO MAKE THIS WORK WITH 2 LAYERS / LENSLETS /ETC
    n_ar_array = np.linspace( (1-tol) * n_ar[0], (1+tol)*n_ar[0], N)
    t_ar_array = np.linspace( (1-tol) * t_ar[0], (1+tol)*t_ar[0], N)
    #print("n_ar_array", n_ar_array)
    
    # LT
    lt_substrate = D['lt_substrate']
    lt_range = np.linspace( (1-tol)*lt_substrate , (1+tol)*lt_substrate, N) 
    print("Range of lt: ", lt_range)

    # NOW WE BEGIN THE MOST UNNEFFICIENT NESTED LOOPS EVER, JUST GO THROUGH ALL PARAMETERS...
    i = 0
    for n_subs in n_substrate_array:
        i=i+1
        print("i: ", i, "N substrate: ", n_subs)  
        for n_ar_test in n_ar_array:
            n_ar_test_array = np.array([n_ar_test])
            
            for t_subs in t_substrate_array:
                
                for t_ar_test in t_ar_array:
                    t_ar_test_array = np.array([t_ar_test]) # this has to be an array
                    
                    for lt_subs in lt_range:
                        for lt_ar_test in lt_range:
                            
                            #print("n_ar_test_array: ", n_ar_test_array)
                            #print("n_ar: ", n_ar)
                            n_array, t_array = SetupLayers(D['opt_type'], n_ar_test_array, N_VACUUM, n_subs, t_subs, t_ar_test_array)
                            lt_array = [lt_ar_test, lt_subs, lt_ar_test, 0]
    
                            #CALCULATE TRANSMISSION
                            Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, t_array, lt_array, freq_array)
                            
                            Ts_array_test = 1 - Refl_s

                            #calculate error
                            error = CalculateError(T_array_data, Ts_array_test)
                            
                            if error < min_error:
                                Ts_array_best = Ts_array_test
                                n_array_best = n_array
                                lt_array_best = lt_array
                                #print("Best lt: ", lt_array_best)
                                t_array_best = t_array
                                min_error = error
                                
    return Ts_array_best, n_array_best, lt_array_best, t_array_best, min_error

######################################

def CalculateError(T_array_data:np.ndarray, Ts_array_test:np.ndarray):
    """Calculates the error between the data transmission and the test, fitted transmission.

    Args:
        T_array_data (np.ndarray): transmission data
        Ts_array_test (np.ndarray): transmission guess

    Returns:
        error [float]: the sum of the absolute error 
    """
    error = np.sum(abs(T_array_data - Ts_array_test))
    return error

####################################################################

def SetupLayers(optics_type:str, n_ar:np.ndarray, n_vacuum:float, n_substrate:float, t_substrate:float, t_ar:np.ndarray):
    """ This function creates the optics layers based on the inputs defined by the user

    Args:
        optics_type (str): [description]
        n_ar (np.ndarray): [description]
        n_vacuum (float): [description]
        n_substrate (float): [description]
        t_substrate (float): [description]
        t_ar (np.ndarray): [description]

    Returns:
        [type]: [description]
    """
    
    if optics_type == "lens":
        n_ar_flipped = np.flipud(n_ar)
        t_ar_flipped = np.flipud(t_ar)
        n_array = np.concatenate( [ [n_vacuum], n_ar, [n_substrate], n_ar_flipped, [n_vacuum] ])
        
        t_array = np.concatenate([ [0], t_ar, [t_substrate], t_ar_flipped])
        
    if optics_type == "coupon":
        n_array = np.concatenate( [ [n_vacuum], n_ar, [n_substrate], [n_vacuum] ])
        t_array = np.concatenate( [ [0], t_ar, [t_substrate] ])
        
    if optics_type == "lenslet":
        n_array = np.concatenate([ [n_vacuum], n_ar, [n_substrate]])
        t_array = np.concatenate( [ [0], t_ar] )   

    return n_array, t_array




#################################
def SelectFrequencyBand(freq_min:float, freq_max:float, freq_array_raw, T_array_raw):
    """Select frequency band to be analyzed by trimming extremities

    Args:
        freq_min (float): if <0, use min(freq_array_raw)
        freq_max (float): if <0, use max(freq_array_raw)
        freq_array_raw ([type]): [description]
        T_array_raw ([type]): [description]

    Returns:
        [type]: [description]
    """
    
    if freq_min < 0:
        freq_min = min(freq_array_raw)
    if freq_max < 0:
        freq_max = max(freq_array_raw)
    
    #find indexes where min and max freq in band appear
    freq = freq_array_raw[0]
    f1 = 0
    while freq < freq_min:
        f1 = f1 + 1
        freq = freq_array_raw[f1]
      
    f2 = f1 - 1 # start where we stopped
    while freq < freq_max:
        f2 = f2 + 1
        freq = freq_array_raw[f2]
    
    freq_array = freq_array_raw[f1:f2]
    T_array = T_array_raw[f1:f2]
    return freq_array, T_array


if __name__ == "__main__":
    main()