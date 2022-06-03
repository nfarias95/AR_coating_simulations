# this script reads the transmission data of multiple files and plots them together
# It also calculates the average transmission in the desired bands
from transmission_data_analysis import *
from useful_functions import plot_95_per_cent_line

import numpy as np
from matplotlib import pyplot as plt
#import matplotlib

PLOT_95_LINE = False

def main():
    
    # set plot font
    FONT_SIZE = 20
    LEGEND_SIZE = 17
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams.update({'font.size': FONT_SIZE})
    plt.rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize

    fig = plt.figure(figsize=(8,6))
    um = "\u03BCm"
    
    # READ FILE
    # insert file path and name here
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/February2022 Angle of Incidence/"
    f_name_array = ["date02072022_transmission_2laerdielectric_theta_var_0.csv", 
                    "date02072022_transmission_2laerdielectric_theta_var_15.csv",
                    "date02072022_transmission_2laerdielectric_theta_var_30.csv",
                    "date02072022_transmission_2laerdielectric_theta_var_45.csv",
                    "date02072022_transmission_2laerdielectric_theta_var_60.csv",
                    "date02072022_transmission_2laerdielectric_theta_var_75.csv",
                    "date02072022_transmission_F3_theta_var_0.csv",
                    "date02072022_transmission_F3_theta_var_15.csv",
                    "date02072022_transmission_F3_theta_var_30.csv",
                    "date02072022_transmission_F3_theta_var_45.csv",
                    "date02072022_transmission_F3_theta_var_60.csv",
                    "date02072022_transmission_F3_theta_var_75.csv"]
    label_array = ["D0", "D15", "D30", "D45", "D60", "D75", 
                   "M0", "M15", "M30", "M45", "M60", "M75"]
    
    angles = [0, 15, 30, 45, 60, 75, 90]

    n_files = len(f_name_array) # number of files being analyzed
    print("Number of files: ", n_files)
    
    bands = np.array([ [77, 224],[77, 224], [77, 224], [77, 224], [77, 224], [77, 224], [77, 224], [77, 224], [77, 224], [77, 224], [77, 224], [77, 224] ] )
    #print("Bands: ")
    #print(bands)
    # initialize array to store transmission data
    trans_band_array = np.zeros(n_files)
    
    # Create an array with different marker types and colors for clarity of plots
    marker_type = ["x", ".", "o", "d", "v", "*", ">"]
    colors = 'r', 'g', 'b', 'o', 'y', 'p'
    
    # Column organization in csv file
    x_column = 0 ; y_column = 1;  z_column = -1
    
    # to calculate the change in transmission due to angle of incidence
    dielectric_variation = np.zeros(int(n_files/2))
    metamaterial_variation = np.zeros(int(n_files/2))
    
    
    plt.figure(1)
    i = -1    
    for f_name in f_name_array:
        i=i+1
        #print("\n\n i: %d, File name: %s" %(  i, f_name))
        f_loc = f_path + f_name 
        
        # GET TRANSMISSION RESULTS
        freq_array, trans_array, var_array, title = read_csv_file(f_loc, x_column, y_column, z_column)
        
        # Calculate transmission in band
        freq_start = bands[i, 0]
        freq_end = bands[i, 1]
        #print("Freq start: ", freq_start)
        #print("Freq end: ", freq_end)
        trans_band_array[i] = calculate_transmission_in_band(freq_start, freq_end, freq_array, trans_array) # average transmission in the band
        
        # calculate the change in transmission
        if i == 0:
            trans_band_array_normal_D = trans_band_array[i] # transmission for angle of incidence = 0, two dielectrics
        elif i == n_files/2:
            trans_band_array_normal_M = trans_band_array[i] # transmission for angle of incidence = 0, metematerial
            
        if i < n_files/2:
            
            dielectric_variation[i] = (trans_band_array_normal_D  - trans_band_array[i])/trans_band_array_normal_D * 100
        else:
            metamaterial_variation[i - int(n_files/2)] = (trans_band_array_normal_M  - trans_band_array[i])/trans_band_array_normal_M *100
                    
        #  ---- PLOT RESULTS  ------
        #marker_i = marker_type[i]
        plt.plot(freq_array, trans_array,  label = label_array[i])
        
    
    # Plot line at 95%
    if PLOT_95_LINE:
        plot_95_per_cent_line()
   
        
    # SETUP PLOT DESIGN
    plt.ylim([0.90, 1.00])
    plt.ylabel("Transmission Power")
    plt.legend()
    plt.xlabel("Frequency [GHz]")
    # plt.show()
    
    
    # PRINT TRANSMISSION RESULTS
    print("\nAverage Transmission:")
    print(trans_band_array)
    # print stuff
    #print("\nDielectric variation: ", dielectric_variation)
    #print("Metamaterial variation: ", metamaterial_variation)
    print("\n Angle | Decrease in Transmission [%] | Decrease in Transmission [%]  | Difference in decrease [%] : ")
    print("       |         (dielectric)         |       (metamaterial)          | (dielectric - metamaterial)")
    for i in range(0, int(n_files/2)):
        print("   %2.d  | %14.2f               |  %14.2f               |  %14.2f          " %(angles[i] , dielectric_variation[i], metamaterial_variation[i], 
                                                    dielectric_variation[i] -metamaterial_variation[i]))
    #print(dielectric_variation - metamaterial_variation )
    print("\n\nEnd of program. \n--------------\n\n")
    

    
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