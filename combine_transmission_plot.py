# this script reads the transmission data of multiple files and plots them together
from transmission_data_analysis import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
def main():
    
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/July2021/"
    f_name_array = ["date07052021_curved_surface_1cyl_s_10_r_30_h_270.csv", "date07052021_corrected_1cyl_s_10_tr_30_th_270.csv"]
    label_array = ["Curved", "Flat"]
    
    f_name_array = ["date07052021_curved_surface_1cyl_s_10_r_30_h_270.csv", "date07122021_curved_surface_solid_dielectric_1layer.csv"]
    label_array = ["Metamaterial", "Solid Dielectric"]
    
    marker_type = ["x", ".", ""]
    n_files = len(f_name_array)
    
    print("Number of files: ", n_files)
    
    x_column = 0
    y_column = 1
    z_column = -1
    
    plt.figure(1)
    i = -1
    for f_name in f_name_array:
        i=i+1
        print("File name: ", f_name)
        f_loc = f_path + f_name 
        
        freq_array, trans_array, var_array, title = read_csv_file(f_loc, x_column, y_column, z_column)
        marker_i = marker_type[i]
        plt.plot(freq_array, trans_array, marker = marker_i, label = label_array[i])
        
        
 
    #plt.ylim([0.80, 1.00])
    plt.ylabel("Transmission", fontname='Arial', size='13')
    plt.legend()
    plt.xlabel("Frequency [GHz]", fontname='Arial', size='13')
    plt.show()
    print("End of program.")
    
if __name__ == "__main__":
    main()