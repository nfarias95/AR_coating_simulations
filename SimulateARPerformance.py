# This code simulates the performance of anti-reflection coating

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import CharlesHill_PB2b_AR_optimize as Charlie

from useful_functions import write_data_to_csv_files, shade_bands

SAVE_DATA = True
output_folder = "C:/Users/nicol/Documents/00Research/Data/ThermalSpray/sims/"
output_filename =  "fillite_thicker_by_2mil.csv" # "thick_fillite_thermal_spray.csv" #"high_fillite_ep_thermal_spray.csv"

PLOT_BANDS_OF_INTEREST = True

def main():
    print("\nThis program simulates the performance of anti-reflection coatings")

    #CONSTANTS
    c = 2.99e8
    n_vacuum = 1 # index of refraction of vacuum
    n_si = 3.4 # index of refraction of silicon

    # DEFINE DESIRED FREQUENCY BANDS
    center_freq = 100e9 # center frequency [Hz]
    center_wl = c/center_freq # associated wavelength [m]

    freq_min = 60e9 # Hz
    freq_max = 180e9 # maximum frequency to be simulated [Hz]
    df = 1e8 # interval between frequencies to be plotted [Hz]

    # DEFINE AR COATING PARAMETERS
    
    # for single layer
    #n_ar = math.sqrt(n_si) # index of refraction of ar layer
    #t_ar = center_wl/4/n_ar # thickness of ar layer [m]

    # Define array of frequencies to be simulated
    num_steps = int((freq_max - freq_min + df) /df)
    freq_array = np.linspace(freq_min, freq_max, num=num_steps) # array containing the frequencies to be plotted.


    # Define array with material layers
    # n_ar = 1.98 ; d_ar = 0.00046228
    # n_al = 3.351 ; d_al = 0.006420556
    # n_array = [1, n_ar, n_al, 1]
    # d_array = [0, d_ar, d_al, 1]
    # lt_array = [0, 0, 0, 0]
    # label = "[vacuum, ar, alumina, ar, vacuum]"

    ep_fillite = 2.3
    t_fillite = 21.8 - 3 # mils
    epsilon_array = np.array( [1, ep_fillite, 5.65, 9.7, 5.65, ep_fillite, 1] )
    n_array = np.sqrt(epsilon_array)
    d_arr_in = np.array([0, t_fillite, 10.6, 254, 10.6, t_fillite, 1]) /1e3
    d_array = d_arr_in * 25.4/1e3
    lt_array = np.zeros(7)
    label = ["vacuum", "fillite", "mixture", "alumina", "mixture", "fillite", "vacuum"]
    label = "vacuum fillite mixture alumina mixture fillite vacuum"

    print("targets: ")
    print(epsilon_array)
    print(n_array)
    print(d_arr_in)
    print(d_array)

    # Simulate transmission
    print("Simulatin transmission...")
    Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, d_array, lt_array, freq_array)

    #Plot results
    plt.figure(1)
    fig, ax = plt.subplots()
    ax.plot(freq_array/1e9, 1-Refl_p)
    ax.set_ylabel("Transmission")
    ax.set_title("Adapted from Charlie, " + label)
    ax.set_ylim(0.8, 1.01)

    if PLOT_BANDS_OF_INTEREST:
        f1 = 75
        f2 = 105
        clr = 'gray'
        shade_bands(ax,  f1, f2, clr)
        f1 = 125
        f2= 165
        shade_bands(ax,  f1, f2, clr)
        


    if SAVE_DATA:
        write_data_to_csv_files(freq_array, 1-Refl_p, "frequency[hz]", "transmission", output_filename, output_folder)


    plt.show()

    print("The end. \n\n")


if __name__ == "__main__":
    main()
    





