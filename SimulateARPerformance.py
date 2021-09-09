# This code simulates the performance of anti-reflection coating

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import CharlesHill_PB2b_AR_optimize as Charlie

def main():
    print("\nThis program simulates the performance of anti-reflection coatings")

    #CONSTANTS
    c = 2.99e8
    n_vacuum = 1 # index of refraction of vacuum
    n_si = 3.4 # index of refraction of silicon

    # DEFINE DESIRED FREQUENCY BANDS
    center_freq = 100e9 # center frequency [Hz]
    center_wl = c/center_freq # associated wavelength [m]

    freq_min = 10e9 # Hz
    freq_max = 150e9 # maximum frequency to be simulated [Hz]
    df = 1e8 # interal between frequencies to be plotted [Hz]

    # DEFINE AR COATING PARAMETERS
    
    # for single layer
    n_ar = math.sqrt(n_si) # index of refraction of ar layer
    t_ar = center_wl/4/n_ar # thickness of ar layer [m]

    # Define array of frequencies to be simulated
    num_steps = int((freq_max - freq_min + df) /df)
    freq_array = np.linspace(freq_min, freq_max, num=num_steps) # array containing the frequencies to be plotted.


    # Define array with material layers
    n_array = [n_vacuum, n_ar, n_si] # index of refraction of materials
    label = "[vacuum,  ar, si]"
    d_array = [0, t_ar] # thickness of material
    lt_array = np.array([0, 0]) # loss (tangent?)

    # Simulate transmission
    print("Simulatin transmission...")
    Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = Charlie.calc(n_array, d_array, lt_array, freq_array)

    #Plot results
    plt.figure(1)
    plt.plot(freq_array/1e9, 1-Refl_p)
    plt.ylabel("Transmission")
    plt.title("Adapted from Charlie, " + label)
    plt.ylim(0, 1.1)

    plt.show()

    print("The end. \n\n")


if __name__ == "__main__":
    main()
    





