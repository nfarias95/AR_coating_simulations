# this script reads the transmission data of multiple files and plots them together
# It also calculates the average transmission in the desired bands
from transmission_data_analysis import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
def main():
    
    # READ FILE
    # insert file path and name here
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/July2021/"    
    f_name_array = ["FB1.csv", "FB2.csv", "FB3.csv"]
    label_array = ["FB1", "FB2", "FB3"] # corresponding label of each file
    n_files = len(f_name_array) # number of files being analyzed
    print("Number of files: ", n_files)
    
    # FREQUENCY BANDS
    # here we define the desired frequency band of each transmission file. 
    # row: file
    # column 0: left limit of band  [GHz]
    # column 1: right limit of band [GHz]
    # this is used to calculated the average transmission
    bands = np.array([ [34, 99] , [60, 162] , [77, 224] ])
    print("Bands: ")
    print(bands)
    # initialize array to store transmission data
    trans_band_array = np.zeros(n_files)
    
    # Create an array with different marker types for clarity of plots
    marker_type = ["x", ".", ""]
    
    
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
        
        # GET TRANSMISSION RESULTS
        freq_array, trans_array, var_array, title = read_csv_file(f_loc, x_column, y_column, z_column)
        
        # Calculate transmission in band
        freq_start = bands[i, 0]
        freq_end = bands[i, 1]
        print("Freq start: ", freq_start)
        print("Freq end: ", freq_end)
        trans_band_array[i] = calculate_transmission_in_band(freq_start, freq_end, freq_array, trans_array)
        
        # PLOT RESULTS
        marker_i = marker_type[i]
        plt.plot(freq_array, trans_array, marker = marker_i, label = label_array[i])
        
        
    # SETUP PLOT DESIGN
    plt.ylim([0.90, 1.00])
    plt.ylabel("Transmission", fontname='Arial', size='13')
    plt.legend()
    plt.xlabel("Frequency [GHz]", fontname='Arial', size='13')
    plt.show()
    
    
    # PRINT TRANSMISSION RESULTS
    print("Average Transmission:")
    print(trans_band_array)
    
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
    
    while freq <= freq_end:
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