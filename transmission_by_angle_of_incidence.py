# this script reads the transmission data of multiple files and plots them together
# It also calculates the average transmission in the desired bands
from transmission_data_analysis import *
from useful_functions import plot_95_per_cent_line

from scipy import integrate
import numpy as np
from matplotlib import pyplot as plt
#import matplotlib

degree_sign = u'\N{DEGREE SIGN}'
PLOT_95_LINE = False
pi = np.pi


# set plot font, other formatting things
FONT_SIZE = 20
LEGEND_SIZE = 17
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': FONT_SIZE})
plt.rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize
#fig = plt.figure(figsize=(8,6))
um = "\u03BCm"

# desired bandwidths:
F1 =  [34,99]
F2 = [60, 162]
F3 = [77, 224]

def main():
    
    # Create an array with different marker types and colors for clarity of plots
    marker_type = ["x", ".", "o", "d", "v", "*", ">"]
    colors = 'r', 'g', 'b', 'o', 'y', 'p'
    # Column organization in csv file
    x_column = 0 ; y_column = 1;  z_column = -1
    
    #  ---- READ FILES -----
    #insert file path and name here
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/February2022 Angle of Incidence/varying_height_finer/"
    f_name_roots = ["date02142022_transmission_F3_theta_var_", "date02142022_transmission_F3_h_plus5_theta_var_",
                    "date02142022_transmission_F3_h_plus10_theta_var_", "date02142022_transmission_F3_h_plus20_theta_var_"]
    root_label = ["Nominal", "+5%", "+10%", "+20%"]
    bands = np.array([F3, F3, F3, F3])
        
    # choose angles to be analyzed
    increment = 5
    max_angle = 65
    angles = np.linspace(0, max_angle, num=int(max_angle/increment+1), dtype=int)
    my_title = "Effect of angle of incidence on mean transmission"
    
    
    test_array = np.array([])
    
    ii=-1
    for f_name_root in f_name_roots:
        ii = ii + 1
        
        band = bands[ii]
        print("\nBand: ", band)
        print(root_label[ii])
        
        label_array = angles
        f_name_array, n_files = get_file_names(angles, f_name_root)
        #print("Number of files: ", n_files)
        
        # initialize array containing the average transmission in the desired band
        trans_band_array = np.zeros(n_files)

        i = -1    
        for f_name in f_name_array:
            i=i+1
            #print("\n\n i: %d, File name: %s" %(  i, f_name))
            f_loc = f_path + f_name 
            
            # GET TRANSMISSION RESULTS
            freq_array, trans_array, var_array, title = read_csv_file(f_loc, x_column, y_column, z_column)
            
            # Calculate transmission in band
            freq_start = band[0]
            freq_end = band[1]
            #print("Freq start: ", freq_start)
            #print("Freq end: ", freq_end)
            trans_band_array[i] = calculate_transmission_in_band(freq_start, freq_end, freq_array, trans_array) # average transmission in the band
            
            # calculate the change in transmission
            if i == 0:
                trans_band_array_normal_M = trans_band_array[i] # transmission for angle of incidence = 0
                        
        
        #print("\ntrans_band_array:", trans_band_array)
        # calculate transmission if we use half half of height
        if ii == 0:
            test_array = np.append(test_array, trans_band_array[0:int(n_files/2)])
            #print("Hey")
            #print("test array: ", test_array)
        elif ii==1:
            test_array = np.append(test_array, trans_band_array[int(n_files/2):])
            #print("ho")
            #print("test array: ", test_array)
        # PLOT AVERAGE TRANSMISSION OVER ANGLE OF INCIDENCE
        plt.figure(2, figsize=(8,6))
        plt.plot(angles, trans_band_array, label=root_label[ii], marker=marker_type[ii])

        
        #  ------  INTEGRATE OVER THETA TO GET THE TOTAL TRANSMISSION IN SPHERE ------
        # dS = r^2 * sin(theta) * d(phi) * d(theta)
        # <T> = 1/Si * integral of surface of Transmission(theta, phi) dphi dtheta
        angles_rad = angles * pi/180 # angles in radians
        Tave = nini_transmission_integral(angles_rad, trans_band_array)
        
        # --  let's integrate using numpy and see how it goes ---
        x = angles_rad
        y = np.multiply(trans_band_array , np.sin(angles_rad) )
        Tave_numpy = np.trapz(y, x) * (np.pi/2 / angles_rad[-1])
        
        # print results
        #print("\nHole depth: ", root_label[ii])
        print("Transmission at normal incidence: %.2f" %(trans_band_array[0]*100) )
        print("Transmission simple average from 0 to %.f degrees: %.2f" %(max_angle, np.mean(trans_band_array*100)))
        #print("Transmission integrated average: %.2f" %(Tave*100))
        #print("Transmission trapezoidal integrated average: %.2f" %(Tave_numpy*100))
    
    
    
    #print("angles: ", angles)
    #print("test array; ", test_array)
    #plt.plot(angles, test_array, label="combined")
    #print("Mean of combined array: ", np.mean(test_array))
    
    
    # format plot
    #plt.figure(2, figsize=(8,6))
    plt.xlabel("Angle of Incidence [Deg]")
    plt.ylabel("Mean Transmittance")
    current_values=plt.gca().get_yticks()
    plt.gca().set_yticklabels(['{:,.2f}'.format(y) for y in current_values])
    
    #plt.title(my_title)
    plt.legend()
    
    
    
    plt.show()
    print("\n\nEnd of program. \n--------------\n\n")
    
def get_file_names(angles, f_name_root):
    f_name_array = []
    # get all file names
    for i in range(0, len(angles)):
        angle = angles[i]
        f_name = f_name_root + str(angle) + ".csv"
        f_name_array.append(f_name)
    
    n_files = len(f_name_array) # number of files being analyzed
    return f_name_array, n_files
    

def nini_transmission_integral(angles_rad, trans_band_array):
    max_angle = angles_rad[-1] # max angle of incidence [radians]
    #print("Max angle:", max_angle*180/np.pi)
    area_fraction = max_angle/(pi/2) # area of the hemisphere explored
    Tsum=0
    for i in range(1, len(angles_rad)):
        d_theta = angles_rad[i] - angles_rad[i-1] # get the angle increment
        #print("dtheta: ", d_theta*180/np.pi)
        #dT = trans_band_array[i] * np.sin(angles_rad[i]) * d_theta / area_fraction
        #dT = d_theta * trans_band_array[i]  * np.sin( angles_rad[i])      
        dT = d_theta * 1  * np.sin( angles_rad[i])   # use this to sanity check the integral      
        Tsum = Tsum + dT
    
    Tave = pi/2 / max_angle * Tsum # divide by the ratio of the total area integrated
    #print("\nTave: ", Tave)
    return Tave
    
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