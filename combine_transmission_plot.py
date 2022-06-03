# this script reads the transmission data of multiple files and plots them together
# It also calculates the average transmission in the desired bands
from transmission_data_analysis import *
import numpy as np
from matplotlib import pyplot as plt
from useful_functions import shade_bands
#import matplotlib


F1 =  [34,99]
F2 = [60, 162]
F3 = [77, 224]

def main():
    
    # set plot font
    FONT_SIZE = 20
    LEGEND_SIZE = 17
    COLOR_BAND_REGION = True
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams.update({'font.size': FONT_SIZE})
    plt.rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize

    #fig = plt.figure(figsize=(8,6.5))

    um = "\u03BCm"
    
    # READ FILE
    # insert file path and name here
    #f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/July2021/"    
    #f_name_array = ["FB1.csv", "FB2.csv", "FB3.csv"]
    #label_array = ["LF-12", "LF-34", "MF-12"] # corresponding label of each file
    #bands = np.array([ F1, F2, F3 ])#), F1, F1, F1 , F1] ) # it doesn't really matter in this case


    f_path = "C:/Users/nicol/Documents/00Research/Data/ThermalSpray/sims/"
    f_name_array = ["ideal_thermal_spray.csv", "high_fillite_ep_thermal_spray.csv", "fillite_thicker_by_2mil.csv", "very_thick_fillite_thermal_spray.csv"]
    label_array = ["ideal", "high n", "thicker by 2 mil", "thicker by 8 mil"]

    # f_path = 'C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/October_2021_redoing_conical_sims/'
    # f_name_array = ["500 μmdate09292021_conical_TR_140_TH_var_s_10.csv",
    #                 "1000 μmdate09292021_conical_TR_140_TH_var_s_10.csv",
    #                 "1500 μmdate09292021_conical_TR_140_TH_var_s_10.csv",
    #                 "2000 μmdate09292021_conical_TR_140_TH_var_s_10.csv",
    #                 "2500 μmdate09292021_conical_TR_140_TH_var_s_10.csv",
    #                 "3000 μmdate09292021_conical_TR_140_TH_var_s_10.csv",
    #                 "K_taper_output_file.csv"]
    # label_array = ["500 " + um, "1000 " + um, "1500 " + um, "2000 " + um, "2500 " + um, "3000 " + um, "KT 2800 " + um]
    
    # FREQUENCY BANDS
    # here we define the desired frequency band of each transmission file. 
    # row: file
    # column 0: left limit of band  [GHz]
    # column 1: right limit of band [GHz]
    # this is used to calculated the average transmission
    bands = [ [75, 105], [125, 165]]
    print("Bands: ")
    print(bands)
    
    n_files = len(f_name_array) # number of files being analyzed
    print("Number of files: ", n_files)
    # initialize array to store transmission data
    trans_band_array = np.zeros(n_files)
    
    # Create an array with different marker types for clarity of plots
    marker_type = ["x", ".", "o", "d", "v", "*", ">"]
    color_array = ['red', 'blue', 'green', 'purple', 'yellow', "violet", "brown"]
    
    
    # Column organization in csv file
    x_column = 0
    y_column = 1
    z_column = -1
    
    plt.figure(1)
    i = -1
    
    for f_name in f_name_array:
        i=i+1
        print("File name: ", f_name)
        f_loc = f_path + f_name 
        
        clr = color_array[i]
        
        # GET TRANSMISSION RESULTS
        freq_array, trans_array, var_array, title = read_csv_file(f_loc, x_column, y_column, z_column)
        
        # Calculate transmission in band
        # freq_start = bands[i, 0]
        # freq_end = bands[i, 1]
        # print("Freq start: ", freq_start)
        # print("Freq end: ", freq_end)
        # trans_band_array[i] = calculate_transmission_in_band(freq_start, freq_end, freq_array, trans_array)
        
        # PLOT RESULTS
        marker_i = marker_type[i]
        plt.plot(freq_array/1e9, trans_array, marker = marker_i, markersize = 1, label = label_array[i], color=clr)
    
        
    # shade the bands region
    ax = plt.gca()
        
    if COLOR_BAND_REGION:
        for band in bands:
            shade_bands(ax,  band[0], band[1], clr)
        
        
    # Plot line at 95%
    xmin = 60
    xmax = 180
    x95 = np.linspace( xmin, xmax, 3)
    y95 = 0.95 * np.ones(3)
    plt.plot(x95, y95, '--', markersize=0.00001, color='k')    
        
    # SETUP PLOT DESIGN
    plt.ylim([0.8, 1.00])
    plt.xlim([xmin, xmax])
    plt.ylabel("Transmittance")
    plt.legend()
    plt.xlabel("Frequency [GHz]")
    plt.show()
    
    
    # PRINT TRANSMISSION RESULTS
    # print("Average Transmission:")
    # print(trans_band_array)
    
    print("\nEnd of program. \n--------------\n\n")
    

    
def calculate_transmission_in_band(freq_start:float, freq_end:float, freq_array:np.ndarray, trans_array:np.ndarray):
    """This function calculate the average transmission in the specified frequency band

    Args:
        freq_start ([type]): [description]
        freq_end ([type]): [description]
        freq_array ([type]): [description]
        trans_array ([type]): [description]
    """
    
    # find index location of bands in freq_array
    freq = 0
    loc = 0
    while(freq < freq_start):
        freq = freq_array[loc]
        loc = loc + 1
    
    loc1 = loc
    
    while freq < freq_end:
        freq = freq_array[loc]
        loc = loc + 1
    
    loc2 = loc
    
    #print("Loc 1:", loc1, "Loc 2: ", loc2)
    # calculate average transmission
    average_trans = np.average(trans_array[loc1:loc2])
    
    #print("Average transmission: ", average_trans)
    return(average_trans)
    
    
    return 0
    
if __name__ == "__main__":
    main()