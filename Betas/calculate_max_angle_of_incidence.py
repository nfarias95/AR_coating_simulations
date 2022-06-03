# Author: Nicole Farias
# Date: July 2021
# Description: This code calculates the transmission through a surface (located in the z plane)
# from the directivity data provided by HFSS
# It assumes that any gain from theta < 90 is reflected power, whereas transmitted power comes from 
# gains with theta > 90


#Importing relevant libraries
import csv
import numpy as np
from matplotlib import pyplot as plt
import math

pi = math.pi
sin = math.sin
cos = math.cos

# Main function
def main():
    print("\n\n Calculate max angle of incidence")
    
    R = 6 #  radius of PIXEL -- assuming spherical cap
    n_si = 3.4
    n_ar1 = 2.5
    n_ar2 = 1.4
    #CHECK ON THIS
    ratio = 0.4 # if I remember correctly, it's R/ext?
    L = R * ratio # extension
    
    # assuming lens isn't spherical cap
    #phi = np.arctan(extension/R)
    
    #theta = pi/2 - phi
        
    # thetas = []
    # phis = np.linspace(0, pi/2, num=30)
    # print(" Phi  theta_r_1, theta_i_1, theta_r_2, theta_i_2, theta_r_3, theta_i_3 ")
    # for phi in phis:
        
    #     theta_r_si, theta_i_si = get_incident_angle_silicon_only(phi, R, L, n_si)
        
    #     theta_r_1, theta_i_1, theta_r_2, theta_i_2, theta_r_3, theta_i_3 = get_incident_angle_two_ARs(phi, R, L, n_si, n_ar1, n_ar2)
        
    #     #print(phi * 180/pi, theta_r * 180/pi, theta_i * 180/pi) 
    #     print("%3.1f   |  %3.1f     %3.1f   |  %3.1f     %3.1f    %3.1f   %3.1f    %3.1f    %3.1f" %(phi*180/pi, 
    #         theta_r_si*180/pi, theta_i_si*180/pi,
    #         theta_r_1*180/pi, theta_i_1*180/pi, theta_r_2*180/pi, theta_i_2*180/pi, theta_r_3*180/pi, theta_i_3*180/pi) )
        
    #     thetas.append(theta_i_3)
    
    """
    thetas_i = np.linspace(0, pi/2, num=10)
    
    for theta_i in thetas_i:
        x = get_antenna_max_size(theta_i, n_si, R, L)
        
        print("x: %.2f [mm] , x/R:  %.4f " %(x, x/R))
    """
    
    for ratio in [0.36, 0.38, 0.39, 0.4, 0.41, 0.42,  0.44]:
        my_L = R*ratio
        calculate_internal_reflection_things(n_si, R, my_L)

    
    print("The end")
    
def degree(theta):
    return theta*180/pi

def rad(theta):
    return theta * pi/180
    
def calculate_internal_reflection_things(n_si:float, R:float, L:float):
    # n_si: index of refraction of silicon
    # R : radius of lenslet
    # L : lenslet extension length
    
    theta_crit = np.arcsin(1/n_si) # critical angle for reflection
    
    thetas = np.linspace(rad(30), rad(80), 500)
    
    for theta in thetas:
        theta_prime = np.arctan(R*sin(theta)/(L+R*cos(theta)))
        theta_i = theta-theta_prime # angle of incidence
        if theta_i > theta_crit:
            print("\n\nL/R : %.2f" %(L/R))
            print("Theta: %.1f , Theta'= %.1f , Theta_i = %.1f (Theta_crit) = %.1f" %(degree(theta), 
                    degree(theta_prime), degree(theta_i), degree(theta_crit)))
            break
    
def get_antenna_max_size(theta_i, n_si, R, L):
    # I am not sure this is working, not accounting for internal reflection
    theta_r = np.arcsin(1/n_si * np.sin(theta_i))
    
    #print("Theta r: %.2f " %(theta_r*180/pi))
    phi = pi/2 - theta_i
    
    psi = pi - phi - theta_r
    
    print("\nThetai: %.2f Theta r: %.2f   phi: %.2f   psi: %.2f" %(theta_i*180/pi, theta_r*180/pi, phi*180/pi, psi*180/pi))
    
    l = np.sin(theta_r) * R / np.sin(psi)
    
    alpha = psi - pi/2
    
    lprime = L * np.tan(alpha)
    
    print("l: %.2f , alpha:  %.2f  , lprime:  %.2f  " %(l, alpha * 180/pi, lprime) )
    
    x = l - lprime
    
    return(x)
    

def get_incident_angle_silicon_only(phi, R, L, n_si):
    P = np.sqrt(R**2+ L**2 + 2*np.cos(phi + np.pi/2))
    # interface between bottom AR and silicon
    theta_r_1 = np.arcsin(L/P * np.sin(phi + np.pi/2))
    print()
    theta_i_1 = np.arcsin(n_si * np.sin(theta_r_1))
    return(theta_r_1, theta_i_1)
    
def get_incident_angle_two_ARs(phi, R, L, n_si, n_ar1, n_ar2):
    
    P = np.sqrt(R**2+ L**2 + 2*np.cos(phi + np.pi/2))
    
    # interface between bottom AR and silicon
    theta_r_1 = np.arcsin(L/P * np.sin(phi + np.pi/2))
    theta_i_1 = np.arcsin(n_si/n_ar1 * np.sin(theta_r_1))
    
    # interface between bottom AR and top AR
    theta_r_2 = theta_i_1 # roughly (assuming R>> AR thickness)
    theta_i_2 = np.arcsin(n_ar1/n_ar2 * np.sin(theta_r_2))
    
    # interface between top AR and vacuum
    theta_r_3 = theta_i_2 # roughly
    theta_i_3 = np.arcsin(n_ar2 * np.sin(theta_r_3))
    
    return(theta_r_1, theta_i_1, theta_r_2, theta_i_2, theta_r_3, theta_i_3)
    

def circ(r):
    # r: radius
    # get the circumference
    return 2*np.pi*r
    
    

if __name__ == "__main__":
    main()