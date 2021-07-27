# Author: Nicole Farias
# Date: July 2021
# Description: This code calculates the transmission through a surface (located in the z plane)
# from the directivity data provided by HFSS for multiple frequencies and plots the transmission over freq.
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
    
    print("\n\n-------------\nWelcome. Let's calculate transmissions based on directivity data!")
    
    #INSERT FILE LOCATION HERE!
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/July2021/"
    f_name = "date07262021_directivity_angle0_silicon_sweep_158_160GHz_box6mm_beam2mm_reformatted.csv"
    f_loc = f_path + f_name
    
    #READ THE HFSS DATA
    freq_array_raw, phi_array_raw, theta_array_raw, dir_array_raw = read_csv_file_freq(f_loc)
    
    
    print("Min|Max freq : ", min(freq_array_raw), max(freq_array_raw))
    print("Min|Max theta : ", min(theta_array_raw), max(theta_array_raw))
    print("Min|Max phi   : ", min(phi_array_raw), max(phi_array_raw))
    
    # Total radiated power (insert here!):
    Prad = 1.0 # [W] total radiated power, in Watts   -- this is not doing much because I actually don't know what Prad is. HFSS let's me 
    #pick the maximum electric field value of the gaussian beam but not the power... 
    
    #Calculate transmission across entire sphere for each different frequency
    freq_array = np.unique(freq_array_raw) # array containing unique frequency values
    nfreqs = len(freq_array)     # number of frequencies
    npoints = int(len(phi_array_raw) / nfreqs) # number of points per frequency
    print("Number of points per frequency: ", npoints)
    print("Freq array: ", freq_array, "\n---------------------")
    
    #Create arrays to store data
    trans_array = np.zeros(nfreqs)
    refl_array = np.zeros(nfreqs)
    
    #LOOP
    
    for counter in range(0, nfreqs):
        freq = freq_array[counter]
        print("\n\nFrequency: ", freq, "GHz")
        
        # get the phis, thetas and directivity values
        if counter == 0: # only need to do this once
            phi_array = phi_array_raw[counter*npoints: (counter+1)*npoints]
            theta_array = theta_array_raw[counter*npoints: (counter+1)*npoints]
        dir_array = dir_array_raw[counter*npoints: (counter+1)*npoints]
        
        #calculate transmission
        trans, refl = calculate_transmission_on_sphere(Prad, npoints, phi_array, theta_array, dir_array)
        trans_array[counter] = trans
        refl_array[counter] = refl
        
    print("Trans array: ", trans_array)
        
    #PLOT RESULTS
    plt.plot(freq_array, trans_array)
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Transmission")
    plt.show()

    print("\nEnd of program. \n------------\n\n")
    

#################################
def calculate_transmission_on_sphere(Prad:float, npoints:int, phi_array:np.ndarray, theta_array:np.ndarray, dir_array:np.ndarray):
    
    #go through each angle and calculate total dir above and below z plane
    Ptrans = 0 # amount of power transmitted
    Prefl = 0 # amount of power reflected 
    
    # Calculate angle increment
    delta_phi = phi_array[1] - phi_array[0]
    print("Delta phi   : ", delta_phi, "degrees")
    points_to_360 = int(360/delta_phi) # how many points until theta goes around 360 degrees 
    delta_theta = theta_array[points_to_360 + 1 ] - theta_array[0] # phi updates everytime theta increases by 360 degrees
    print("Delta theta : ", delta_theta, "degrees")
    
    delta_phi_rad = delta_phi * pi/180
    delta_theta_rad = delta_theta * pi/180
    
    
    # go through all solid angles and find directivity
    # if theta is above 90, light is transmitted. Otherwise it is reflected
    for n in range(0, npoints):
        
        theta = theta_array[n] # angle theta 
        theta_rad = theta * pi/180
        
        d = dir_array[n] # directivity at some angle (theta, phi)
        U = Prad * d/(4*pi) # radiation intensity (W/unit solid angle)
        
        deltaP = U * math.sin(theta_rad) * delta_theta_rad * delta_phi_rad # delta Power (power at that solid angle)
        
        #CALCULATE REFLECTED POWER
        if theta < 90 : # above the z-plane
            Prefl = Prefl + deltaP
        
        # CALCULATE TRANSMITTED POWER
        else:
            Ptrans = Ptrans + deltaP
        
    Ptotal = Ptrans + Prefl
    print("Transmitted Power: %.2f, Reflected power: %.2f, Total power: %.2f [W]"  %( Ptrans, Prefl, Ptotal ))
    
    # Calculate reflection and transmission
    trans = Ptrans/Ptotal
    refl = Prefl/Ptotal
    
    print("Percentages: Transmitted: %.3f    Reflected: %.3f " % (trans, refl) )
    return(trans, refl)
        
        
        
        
    
#################################  
def read_csv_file_freq(f_loc:str):
    """ This function reads a csv file

    Args:
        f_loc ([type]): file location

    Returns:
        phi (degrees), theta (degrees),  dir total 
    """
    
    
    fields = []
    rows = []
    
    with open(f_loc, 'r') as csvfile:
        print("Reading data")
        #creating a csv reader object
        csvreader = csv.reader(csvfile)
        
        #extracting field names through first row
        fields = next(csvreader)
        
        #extracting each data row one by one
        for row in csvreader:
            rows.append(row)
            
        #get total number of rows
        print("\nTotal number of rows: %d " % (csvreader.line_num))
        
    # printing th field names
    print('\nFile header: ' + ', ' .join(field for field in fields))
    
    # printing first 5 rows
    print("\n First 5 rows of data file are: \n")
    counter = -1
    for row in rows[:5]:
       counter = counter + 1
       print(rows[counter])
       #for col in row:
           #print("%10s " % col)
       print('\n')
        
        
    freq_array = np.zeros(csvreader.line_num -1) # array to store phi angle
    phi_array = np.zeros(csvreader.line_num -1) # array to store phi angle
    theta_array = np.zeros(csvreader.line_num -1) # array to store theta angle
    dir_array = np.zeros(csvreader.line_num -1) # array to store directivity data at a specific phi and theta
    
    # put data into each array
    counter = -1
    for row in rows:
        counter = counter + 1
        freq_array[counter]  = float(row[0])
        phi_array[counter]   = float(row[1])
        theta_array[counter] = float(row[2])
        dir_array[counter]   = float(row[3])
    
    
    
    return freq_array, phi_array, theta_array, dir_array
    
    
if __name__ == "__main__":
    main()