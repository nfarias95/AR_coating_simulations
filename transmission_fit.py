#This program does a simple fit of transmission data to get the thickness and index of refraction of anti-reflection coatings
# To do: import parameters from a text file. Right now things are too confusing!


import math
import numpy as np
import cmath
import matplotlib.pyplot as plt
import CharlesHill_PB2b_AR_optimize as Charlie
from transmission_data_analysis import *


def main():

    print("Welcome! This program calculates the index of refraction and thickness of anti-reflection coatings \nby doing a simple fit of transmission data. \nYou must provide the data file and the ranges of interest.")

    # DATA FILE LOCATION
    filepath = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/July2021/"
    filename = "date07052021_corrected_1cyl_s_10_tr_var_th_475.csv"
    f_loc = filepath + filename
    
    #select columns where data is located
    freq_column = 2
    data_column = 3
    var_column = 1 # column containing a variation parameter (ex: radius, height, spacing). For simulations only. Set -1 if no parameters are varying
    var_selection = 0 # select which of the variations we are looking at.
    var_label = "?"
    
    #Read data file
    freq_array_raw, T_array_raw, var_array, title = read_csv_file(f_loc, freq_column,  data_column, var_column)
    #split data if necessary
    if var_column >= 0:
             freq_array_raw, T_array_raw, var_values  = split_data(var_array, freq_array_raw, T_array_raw)
             #select one of the columns
             T_array_raw = T_array_raw[:, var_selection]
             var_label = str(var_values[var_selection])
    else:
        var_values = []

    #CONSTANTS
    #indices of refracion
    n_vacuum = 1 # index of refraction of vacuum
    n_substrate = 3.11 # index of refraction of substrate/lens. Typically known. 3.11 for alumina, 3.4 for silicon
    lt_substrate = 0.00001 # loss tangent of substrate/lens

    #RANGE OF EXPECTED INDEX OF REFRACTION
    n_range = np.linspace(1.2, 3.2, num=40)
    #RANGE OF EXPECTED LOSS TANGENT
    lt_range = np.linspace(0, 1e-6, num=20)
    
    # angle of incidence
    theta0 = 0 

    # thickness of samples:
    t_substrate = 0.25 * 2.54/100 # [m] thickness of alumina witness sample
    
    #t_ar = np.array([0.017, 0.0105]) * 2.54/100 # [m] thickness of AR layers (provided by Oliver)
    t_ar = np.array([475e-6])
    # Initial guess of index of refraction of ar coating
    #n_ar = np.array([math.sqrt(2), math.sqrt(5.2)])
    n_ar = [1.8]

    # SETUP LAYERS
    optics_type = "lenslet" # either "lens" or "lenslet"
    n_array, t_array = SetupLayers(optics_type, n_ar, n_vacuum, n_substrate, t_substrate, t_ar)
               
        
    # Initialize array of loss tangent
    lt_array = np.zeros(len(n_array))

    print("Indices of refraction: \n", n_array)
    print("Thicknesses: \n", t_array)
    
    
    # SELECT FREQUENCY BAND TO BE ANALYZED
    freq_array_raw = freq_array_raw * 1e9 
    freq_min = 20e9 # [Hz]
    freq_max = 200e9 # [Hz]
    freq_array, T_array_data = SelectFrequencyBand(freq_min, freq_max, freq_array_raw, T_array_raw)
    
    # CHECK DATA
    #print("Frequency: \n", freq_array)
    #print("Transmission: \n", T_array_data)
    #plot_results(freq_array, T_array_data, [], "Sim data", [freq_min, freq_max])
    
    print("\n\nn_array: \n", n_array)
    print("t_array: ", t_array)
    print("lt_array: ", lt_array)
    
    #CALCULATE INITIAL GUESS
    #Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, t_array, lt_array, freq_array)


    #CALCULATE BEST FIT
    Ts_array_fit, n_array_fit, lt_array_fit, error = FitDielectricConstant( T_array_data, freq_array, n_array, lt_array, t_array, 
                                                theta0, n_range, lt_range)
    
    
    
    
    #OUTPUT
    print("\n\n=-=-=-=-=-=-=-=-=-= Fit results: =-=-=-=-=-=-=-=-=-=-=")
    print("Error: ", error)
    print("Index of refraction: ", n_array_fit)
    print("Loss tangent: ", lt_array_fit)        
            
    #plot initial guess, data and best fit
    fig2 = plt.figure()
    ax = plt.subplot(111)
    ax.plot(freq_array/(10**9), T_array_data, label = "Data")
    #ax.plot(freq_array/(10**9), Ts_array_baseline, label = "Guess")
    ax.plot(freq_array/(10**9), Ts_array_fit, label = "Fit")
    ax.legend()
    plt.ylabel("Transmission")
    plt.xlabel("Frequency [Ghz]")
    plt.title("TR = " + var_label)
    plt.ylim(0, 1.1)
    plt.show()


    print("\nThe end.\n\n")

###############################################
############# END OF MAIN ####################
###############################################



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
        
    if optics_type == "lenslet":
        n_array = np.concatenate([ [n_vacuum], n_ar, [n_substrate]])
        t_array = np.concatenate( [ [0], t_ar] )   

    return n_array, t_array


def ReadTransmissionData(filepath: str, filename: str):
    """Read the transmission data with format [frequency, transmission, garbage]

    Args:
        filepath (str): path of data file
        filename (str): name of data file

    Returns:
        [type]: [description]
    """
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


def SelectFrequencyBand(freq_min:float, freq_max:float, freq_array_raw, T_array_raw):
    """Select frequency band to be analyzed by trimming extremities

    Args:
        freq_min (float): [description]
        freq_max (float): [description]
        freq_array_raw ([type]): [description]
        T_array_raw ([type]): [description]

    Returns:
        [type]: [description]
    """
    
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


def FitDielectricConstant( T_array_data, freq_array, n_array_0, lt_array_0, t_array, theta0:float, n_range, lt_range):
    """ Calculate a simple best fit of the dielectric constant
    
    Args:
        T_array_data ([type]): array containing transmission data
        freq_array ([type]): array containing frequency data
        n_array_0 ([type]): initial guess of indices of refraction
        lt_array_0 ([type]): initial guess of loss tangent
        d_array ([type]): array containing thickness of materials
        theta0 (float): angle of incidence
        n_range ([type]): range of acceptable indices of refraction
        lt_range ([type]): range of acceptable loss tangents

    Returns:
        [type]: [description]
    """
    #initialize error value
    min_error = 1e5 # start with a large error 
    
    #initialize arrays for best fit results
    Ts_array_best = np.zeros(len(T_array_data))
    n_array_best = np.ones(len(n_array_0))
    lt_array_best = np.zeros(len(lt_array_0))
    
    #initialize index of refraction and loss tangent arrays
    n_array = n_array_0
    lt_array = lt_array_0
    
    #print("n_array", n_array)
    #print("lt_array_best", lt_array_best)
    
    # optimize for each layer
    for layer in range(1, len(n_array) - 1): # the first and last indeces are for vacuum, they don't count
        #print("Layer: ", layer)
        
        #try different loss tangents
        for lt in lt_range:
            #iterate over loss tangent
            #print("\nlt :", lt)
            lt_array[layer - 1] = lt #index is off by one
            
            for n in n_range:
                #iterate over index of refraction
                #print("n: ", n)
                n_array[layer] = n
                
                #CALCULATE TRANSMISSION
                Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, t_array, lt_array, freq_array)
                
                Ts_array_test = 1 - Refl_s
                
                #calculate error
                error = CalculateError(T_array_data, Ts_array_test)
                #print("Minimum Error: %.4f Error: %.4f" %(min_error, error))
                
                #keep the setting with smallest error
                if error < min_error:
                    #print("Better!!!")
                    Ts_array_best = Ts_array_test
                    n_array_best[layer] = n
                    lt_array_best[layer-1] = lt
                    min_error = error
                #print("N array best: ", n_array_best)
                #print("Lt array best: ", lt_array_best)

            
    return Ts_array_best, n_array_best, lt_array_best, min_error


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




if __name__ == "__main__":
    main()

