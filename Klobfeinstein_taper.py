# This code calculates the Klobfeinstein taper

# main equation: ln(Z0(z)) = ln(sqrt(Zs*Zl)) + ln(sqrt(Zl/Zs)) * Gamma_m * (Zl + Zs)/(Zl-Zs) * [A^2 * phi(2x/l, A) + 
#   U(z-1/2*l) + U(z+1/2*l) - 1]

# Source: engineering libre texts, Tapered Matching Transformers

# To do: fix phi function. It is not giving me what I want.

import cmath
import numpy as np
from math import acosh, cosh,  log, sqrt, cos # log is ln?
from scipy.constants import mu_0, epsilon_0, c, pi
from matplotlib import pyplot as plt
import CharlesHill_PB2b_AR_optimize as Charlie
from useful_functions import *
import csv

def main():
    plt.rcParams.update({"font.family": "Times New Roman"})
    um = "\u03BCm"
    
    print("\n-----------------------------------\n\nWelcome! \nLet's calculate the Klobfeinstein taper!")
    
    Zsilicon = calculate_Z_from_n(3.4)
    Zvacuum = calculate_Z_from_n(1)
    print("Impedance of silicon: ", int(Zsilicon))
    print("Impedance of vacuum: ", int(Zvacuum))
    
    
    # CONSTANTS
    npoints = 50 # number of points in impedance calculation
    Zl = Zvacuum
    Zs = Zsilicon
    Gamma_m = 0.05
    epsilon_e = 3.4**2
    freq_cut = 16e9 # to determine the taper lenght
    min_freq = 50e9 # Hz (to plot)
    max_freq = 250e9 # Hz (to plot)
    n_freq_points = 100 # number of points to be plotted in frequency response
    
    # Calculate A
    rho_0 = calculate_rho0(Zl, Zs)
    print("rho_0: ", rho_0)
    A = calculateA(rho_0, Gamma_m)
    print("A: ", A)
      
    
    #find the length of the microstrip:
    #from libre text	length = A * lambda(min freq in band) /(2*pi*sqrt(epsilon_e))
    lambda_max = c/freq_cut #in meters
    l = A * lambda_max / (2*pi*sqrt(epsilon_e) )
    print("Lenght of microstrip: %.2f microns" % (l*1e6) )
    
    z_array = np.linspace(-l/2, l/2, npoints)
    #z_array = np.linspace(0, l, npoints)
    

    # plot phi to see if things are going ok. looks ok
    #plot_phi()
    
    #initialize data arrays
    Z0_array = np.zeros(npoints)
    
    # CALCULATE THE IMPEDANCE AT EACH POINT
    for i in range(0, npoints):
       z = z_array[i]
       #print("\n z: ", z*1e6)
       
       Z0 = calculate_Z0(z, Zs, Zl, rho_0, l, A)

       #print("Z0: ", Z0)
       Z0_array[i] = Z0
    
    # calculate index of refraction at each point
    n_array = calculate_n_from_Z(Z0_array)
    
    plt.figure(1)
    plt.scatter(z_array, Z0_array)
    plt.xlabel("z (microns)")
    plt.ylabel("Z0")
    plt.grid()
    
    plt.figure(2)
    plt.scatter(z_array[0:-1] + l/2, n_array[0:-1])
    plt.xlabel("z " + "um")
    plt.ylabel("n")
    plt.grid()
    
    
    # NOW CALCULATE THE FREQUENCY RESPONSE
    
    transmission_array, freq_array = calculate_transmission_from_Z(Z0_array[0:-1], l, min_freq, max_freq, n_freq_points) # ignore the last point
    plt.figure(4)
    plt.plot(freq_array/1e9, transmission_array)
    plt.ylabel("Transmission")
    plt.title("Transmission of Klobfeinstein taper")
    plt.xlabel("Frequency [GHz]")
    
    # write data to csv file
    folder = 'C:/Users/nicol/Documents/00Research/PythonCode/AR_coating_simulations/output_files/'
    fname = 'K_taper_output_file.csv'
    write_data_to_csv_files(freq_array / 1e9, transmission_array, "freq [GHz]", "transmission", fname, folder)
    
    
    # OPTIMIZE TRANSMISSION IF WE WANT TO TRY THAT
    n1 = 1
    n2 = 3.4
    freq_band = [77e9, 224e9]
    optimize_transmission_in_band(n1, n2, freq_band, Gamma_m)
    
    
    
    #freq_array = np.linspace(min_freq, max_freq, npoints)
    
    # Calculate Beta
    # from excel spreadsheet from https://www.microwaves101.com/encyclopedias/klopfenstein-taper, beta = A/l
    #beta = A/l # rad/cm
    #lambda0 = 2*pi/beta   this kinda weird because they are calculating frequency as opposed to setting a fixed frequency
    
    # let~s try this:
    #lambda_array = c/freq_array
    #beta_array = 2*pi/lambda_array
    
    #reflection_array = np.zeros(len(beta_array))
    
    # from https://www.researchgate.net/publication/325569490_Design_of_a_Klopfenstein_Taper_for_Impedance_Matching_of_a_High-Speed_Photodetector
    # reflection = reflection_0 * sqrt( (beta*l)^2 - (A)^2)/cosh(A)
    # assuming reflection_0 is our rho
    # i=0
    # for beta in beta_array:
    #     refl = rho_0 * cos(sqrt( abs( (beta*l)**2 - (A)**2) ) /cosh(A) )
    #     reflection_array[i] = refl 
        
    #     i=i+1   

    # plt.figure(3)
    # plt.plot(freq_array/1e9, reflection_array)
    # plt.xlabel("Frequency [GHz]")
    # plt.ylabel("Reflection")
    
    plt.show()
    


    print("The end.")
    
    
    
##########################################################
##########################################################
##########################################################

def optimize_transmission_in_band(n1, n2, freq_band, Gamma_m):
    print("Optimizing transmission... \n")
    
    npoints = 50 # number of "layers"
    n_freq_points = 100
    
    freq_start = freq_band[0] * 0.1 # where we start scanning
    
    # calculate the impedances
    Z1 = calculate_Z_from_n(n1)
    Z2 = calculate_Z_from_n(n2)
    
    # calculate epsilon equivalent
    epsilon_e = n2**2
    
    # Calculate A
    rho_0 = calculate_rho0(Z1, Z2)
    print("rho_0: ", rho_0)
    A = calculateA(rho_0, Gamma_m)
    print("A: ", A)
    
    # start going over different start frequencies
    freq = freq_start
    max_ave_trans = 0 # initial value of average transmission
    
    print("\n Starting loop \n")
    
    while freq <= freq_band[0]:
        #find the length of the microstrip:
        #from libre text	length = A * lambda(min freq in band) /(2*pi*sqrt(epsilon_e))
        lambda_max = c/freq #in meters
        l = A * lambda_max / (2*pi*sqrt(epsilon_e) )
    
        z_array = np.linspace(-l/2, l/2, npoints)
        Z0_array = np.zeros(npoints)

        
        
        # CALCULATE THE IMPEDANCE AT EACH POINT. this might be ok outside of loop
        for i in range(0, npoints):
            z = z_array[i]
            Z0 = calculate_Z0(z, Z1, Z2, rho_0, l, A)
            Z0_array[i] = Z0
    
        # calculate transmission
        transmission_array, freq_array = calculate_transmission_from_Z(Z0_array[0:-1], l, freq_band[0], freq_band[1], n_freq_points) # ignore the last point
        ave_trans = np.mean(transmission_array)
        #print("Freq start: ", freq/1e9, "T: ", ave_trans)
        
        if ave_trans > max_ave_trans and l*1e6 < 3000:
            print("Hey! Freq: ", freq/1e9)
            print("l : ", l*1e6)
            max_ave_trans = ave_trans
            best_freq_start = freq
                #Plot results
        
            plt.figure(5)
            plt.plot(freq_array/1e9, transmission_array)
            plt.ylabel("Transmission")
            plt.title("Transmission of Klobfeinstein taper")
            plt.xlabel("Frequency [GHz]")
    
        # update freq
        freq = freq + 0.5e9 # one gigaheartz at a time
        #print("freq: ", freq)
    
    print("Best freq start: ", best_freq_start/1e9)
    
    print("Done Optimizing")
    return(0) # return transmission
    





def calculate_transmission_from_Z(Z_array, l, freq_min, freq_max, n_freq_points):
    
    # initialize frequency array
    freq_array = np.linspace(freq_min, freq_max, n_freq_points) # in Hertz
       
    # get thickness of each layer - uniform "mesh"
    npoints = len(Z_array) # number of "layers"
    t = l/npoints
    t_taper_array = t * np.ones(npoints) # thickness of layers
    t_array = np.append(np.append(0, t_taper_array) , 0) # first and last indices don't matter
    
    # get n of taper from Z
    #print("Zs:", Z_array)
    n_taper_array = calculate_n_from_Z(Z_array)
    n_array = np.append( np.append(3.4, n_taper_array), 1)
    
    
    #print("L: ", l*1e6, " microns")
    #print("Ns:", n_array)
    #print("Ts:", t_array)
    
    # initialize lt array
    lt_array = np.zeros(len(n_array))
    
    # now calculate transmission
    Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, t_array, lt_array, freq_array)
    
    
    return (1-Refl_p, freq_array)




def calculate_n_from_Z(Z):
    # get index of refraction from impedance
    Z0 = np.sqrt(mu_0/epsilon_0)
    n = Z0/Z
    return n



def calculate_Z_from_n(n):
    # get impedance from index of refraction
    Z0 = np.sqrt(mu_0/epsilon_0)
    Z = Z0/n # for non-magnetic media
    return(Z)

def calculate_rho0(Zl, Zs):
    rho_0 = 1/2 * log(Zl/Zs)
    return(rho_0)

    
def calculate_Z0(z, Zs, Zl, rho_0, l, A):
    
    phi = calculate_phi(2*z/l, A) 
       
    lnZ0 = 1/2 * log(Zs*Zl) + rho_0/cosh(A) * ( A**2 * phi + U(z - l/2) + U(z + l/2) ) #  from original paper
       
    #lnZ0 = log(sqrt(Zl*Zs)) + log(sqrt(Zl/Zs)) * Gamma_m * (Zl + Zs)/(Zl-Zs) * ( A**2 * phi + U(z - 1/2*l) + U(z + 1/2*l) -1) # from libre text, not accurate I don't know why
        
    Z0 = np.exp(lnZ0)
    
    #print("ln(Z0): %.2f , Z0: %.2f" %(lnZ0, Z0 ))
    
    return Z0
    
################        
# the function
def U(x):
    if x < 0:
        return(0)
    else:
        return(1)

#############################
def calculate_phi(w, A):
    # the plot phi function shows that this seems to be working well.
    phi = 0 # initialize phi
    k_max = 20 # 20 is sufficient according to website
    
    for k in range(0, k_max):
        # calculate the coefficients
        if k==0:
            ak = 1
            bk = 1/2 * w
            
        else:
            ak = A**2 /( 4 *k * (k+1) ) * akm1
            bk = (1/2 * w * (1-w**2)**k + 2*k*bkm1) / ( 2*k + 1)
        
        # calculate the sum
        phi = phi + ak*bk
        
        # update (k-1) terms
        akm1 = ak # a(k-1)
        bkm1 = bk
    
    #print output?
    #print("Phi: " , phi, "w: ", w)
       
    
    return phi 

def plot_phi():
    A_array = np.array([3.0, 3.5, 5, 5.6])
    w_array = np.linspace(0, 1, 20)
    phi_array = np.zeros(20)
    
    plt.figure(2)
    for A in A_array:
        i = 0
        for w in w_array:
            
            phi = calculate_phi(w, A)
            phi_array[i] = phi
            
            # update counter
            i = i + 1
        
        plt.plot(w_array, phi_array, label="A = " + str(A))
        
    plt.grid()
    plt.xlabel("w")
    plt.ylabel("phi")
    plt.legend()
    
def calculateA( rho_0, Gamma_m):
    A = acosh( rho_0/ Gamma_m)
    return A
    




if __name__ == "__main__":
    main()  