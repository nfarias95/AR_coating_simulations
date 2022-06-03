#!/usr/bin/env python
# coding: utf-8

# Author: Charles Hill (chill90)
# small alterations in the main function by Nicole Farias

import numpy as np
import scipy.optimize as opt
import scipy.constants as ct
import matplotlib.pyplot as plt
import seaborn as sns
import math

def main():
    sns.set(style='whitegrid', font='serif', font_scale=2.5)
    lw = 4
    
    # Constants
    n_vacuum = 1  # index of refraction of vacuum
    n_alumina = 3.11  # index of refraction of alumina. --- > i have to double check this value
    n_silicon = 3.4 # index of refraction of Silicon
    
    center_freq = 100e9
    wl = 3e8/center_freq    
    
    n_metamaterial = math.sqrt(n_silicon) # desired index of refraction of metamaterial layer
    thick_metamaterial = wl/4/n_metamaterial # 675e-6 # thickness of metamaterial [m]
    lenslet_radius = 1.6e-3 # 

    n_arr = [n_vacuum, n_metamaterial, n_silicon]
    d_arr = [0 , thick_metamaterial ]
    lt_arr = np.array([0 , 0])
    
    # FREQUENCY PARAMETERS
    freq_min = 10e9 # Hz
    freq_max = 150e9 # Hz
    df = 1e8  # interval between frequencies to be plotted [Hz]
    num_steps = int((freq_max - freq_min + df)/df)
    freq_arr = np.linspace(freq_min, freq_max, num =num_steps) # array containing the frequencies to be plotted
    
    Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s = calc(n_arr, d_arr, lt_arr, freq_arr)
    
    
    #plot results
    plt.figure(1)    
    plt.plot(freq_arr/1e9, 1 - Refl_p)
    plt.ylabel("Transmission")
    plt.xlabel("Frequency [Ghz]")
    plt.title("Adapted from Charlie")
    plt.ylim(0, 1.1)
    
    plt.show()

# Calculate boundary conditions using the Hou formalism
def calc(n_arr, d_arr, lt_arr, freq_arr, inc_ang=0):
    # Number of interfaces
    num_inters = len(n_arr) - 1
    # Calculate the angles of propogation in each layer
    theta = [inc_ang]
    # Use Snell's law at each interface
    for i in range(num_inters):
        new_ang = np.arcsin((n_arr[i] / n_arr[i + 1])*np.sin(theta[i]))
        theta.append(new_ang)
    theta = np.array(theta)
    # Calculate the reflection coefficients
    rs = np.array([((n_arr[i] * np.cos(theta[i]) - n_arr[i + 1] * np.cos(theta[i + 1])) /
                    (n_arr[i] * np.cos(theta[i]) + n_arr[i + 1]*np.cos(theta[i + 1])))
                    for i in range(num_inters)])
    rp = np.array([((n_arr[i] * (1. / np.cos(theta[i])) - n_arr[i + 1] * (1. / np.cos(theta[i + 1]))) /
                   (n_arr[i] * (1. / np.cos(theta[i])) + n_arr[i + 1]*(1. / np.cos(theta[i + 1]))))
                   for i in range(num_inters)])
    # Calculate transmission coefficients
    ts = np.array([1. + rs[i] for i in range(0, num_inters)])
    tp = np.array([1. + rp[i] for i in range(0, num_inters)])
        
    # Calculate reflection and transmission for each frequency
    Ms = np.matrix([[1., rs[0]], [rs[0], 1.]]) * (1. / ts[0]) #Identity matrix for the s polarization
    Mp = np.matrix([[1., rp[0]], [rp[0], 1.]]) * (1. / tp[0]) #Identity matrix for the p polarization
    for i in range(1, num_inters):
        # Calculate the exponential propogation factors
        pp = (2. * ct.pi * n_arr[i] * d_arr[i] * (
            1. / np.cos(theta[i])) * np.array(freq_arr) / ct.c) * (0.5 * lt_arr[i] + 1.0j)
        pm = (2. * ct.pi * n_arr[i] * d_arr[i] * (
            1. / np.cos(theta[i]))*np.array(freq_arr) / ct.c) * (-0.5 * lt_arr[i] - 1.0j)
        # Calculate new propogation matrices
        Ms = (Ms * np.matrix([[np.exp(pp), 0],[0, np.exp(pm)]]) *
              np.matrix([[1., rs[i]], [rs[i], 1.]])*(1. / ts[i]))
        Mp = (Mp * np.matrix([[np.exp(pp), 0],[0, np.exp(pm)]]) *
              np.matrix([[1., rp[i]], [rp[i], 1.]])*(1. / tp[i]))
            
    # Frequencies
    Freq = freq_arr
        
    # Calculate transmitted power
    Tran_s = np.array(list(map(lambda x: abs(x)**2, 1. / (Ms.item(0, 0)))))
    Tran_p = np.array(list(map(lambda x: abs(x)**2, 1. / (Mp.item(0, 0)))))

    # Calculate reflected power
    Refl_s = np.array(list(map(lambda x: abs(x)**2, Ms.item(1, 0) / Ms.item(0, 0))))
    Refl_p = np.array(list(map(lambda x: abs(x)**2, Mp.item(1, 0) / Mp.item(0, 0))))

    # Calculate absorbed power
    Abso_s = 1. - Tran_s - Refl_s
    Abso_p = 1. - Tran_p - Refl_p
        
    return Freq, Tran_p, Tran_s, Refl_p, Refl_s, Abso_p, Abso_s

def pb2b_trans(n_arr, d_arr, lt_arr, freq_arr, inc_ang=0):
    # Calculate transmission
    out = calc(n_arr, d_arr, lt_arr, freq_arr, inc_ang=0)
    freq = out[0]
    refl_p = out[3]
    refl_s = out[4]
    tran_avg = (refl_p + refl_s) / 2
    # Interpolate the PB2b bands to the given frequencies
    tran_90 = np.interp(freq_90, freq * 1e-9, tran_avg, left=0, right=0)
    tran_150 = np.interp(freq_150, freq * 1e-9, tran_avg, left=0, right=0)
    band_tran_90 = np.trapz(tran_90 * band_90, freq_90) / np.trapz(band_90, freq_90)
    band_tran_150 = np.trapz(tran_150 * band_150, freq_150) / np.trapz(band_150, freq_150)
    return (band_tran_90 + band_tran_150) / 2

def thick_opt_1(p, freq, inc_ang=0):
    n = [1., p[0], n_sapp, p[0], 1.]
    d = [0., p[1], d_sapp, p[1], 0.]
    lt = np.zeros(len(n))
    ret = pb2b_trans(n, d, lt, freq, inc_ang=0)
    return ret

def thick_opt_2(p, freq, inc_ang=0):
    n = [1., p[0], p[1], n_sapp, p[1], p[0], 1.]
    d = [0., p[2], p[3], d_sapp, p[3], p[2], 0.]
    lt = np.zeros(len(n))
    ret = pb2b_trans(n, d, lt, freq, inc_ang=0)
    return ret

def thick_opt_3(p, freq, inc_ang=0):
    n = [1., p[0], p[1], p[2], n_sapp, p[2], p[1], p[0], 1.]
    d = [0., p[3], p[4], p[5], d_sapp, p[5], p[4], p[3], 0.]
    lt = np.zeros(len(n))
    ret = pb2b_trans(n, d, lt, freq, inc_ang=0)
    return ret

    
if __name__ == "__main__":
    main()






