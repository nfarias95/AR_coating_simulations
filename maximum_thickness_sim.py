# This code is an attempt to analyze the effect of thick AR coatings
# Nicole Farias
# Last updated October 2021

import math
from math import sin, cos, asin, acos, tan
import numpy as np
from matplotlib import pyplot as plt
pi = math.pi

def main():
    print("Let's see how thick AR coatings can be")
    
    # CONSTANTS
    n1 = 3.4
    n2 = math.sqrt(n1)
    R = 5.5 # [mm] radius of lenslet
    t = 0.5 # [mm] thickness of single-layer AR coating
    my_title = " n1 = " + str(n1) + ", n2 = " + str(int(n2*100)/100) + ", R = " + str(R) + " mm and t = " + str(t) + " mm"
      
    # make a loop
    npoints = 25
    r0_array = np.linspace(0.1, 4.2, npoints) 
    r_lens_array = np.zeros(npoints) # initialize array
    r_ar_array = np.zeros(npoints) # initialize array
    difference_array = np.zeros(npoints)
    print("r0 array: ", r0_array)
    
    for i in range(0,npoints):
        r0 = r0_array[i]
        
        r_lens = ray_path_lens_only(r0, R+t, n1) # just the lenslet
        r_ar = ray_path_single_ar(r0, R, t, n1, n2) # lenslet with AR machined into it
        
        
        difference_array[i] = abs( (r_lens-r_ar) / r_lens) * 100
        r_lens_array[i] = r_lens
        r_ar_array[i] = r_ar
        
    
    # plot results
    plt.figure(1)
    plt.plot(r0_array, r_lens_array, label="without AR")
    plt.plot(r0_array, r_ar_array, label="with AR")
    plt.xlabel("r0 [mm]")
    plt.ylabel("r [mm]")    
    plt.title(my_title)
    #plt.title("Lenslet Radius = 6mm, AR thickness = 3 mm")
    plt.legend()
    
    plt.figure(2)
    plt.plot(r0_array/(R+t) , difference_array)
    plt.xlabel("r0/(R)")
    plt.ylabel("Percent difference")
    plt.title("Lenslet Radius = 6mm, AR thickness = 3 mm")
    
    
    #plot projection
    plt.figure(4)
    plot_projection(R, r0_array, 'g.')
    plot_projection(R, r_lens_array, 'r.')
    plot_projection(R, r_ar_array, 'b.')
    plt.title(my_title)
    
    plt.show()
    
    print("The end.")
    
    
#################################### 

def plot_projection(R, r_array, color):
    
    n_angles = 100
    theta_array = np.linspace(0, 2*pi - 2*pi/n_angles, n_angles)
    
    for theta in theta_array:
        plt.polar(theta, R, 'g.')
        
        for r in r_array:
            plt.polar(theta, r, color)

    
    return 0
    

def degrees(angle_in_radians):
    return angle_in_radians * 180/pi  
  
def ray_path_lens_only(r0:float, R:float, n1:float):
    """[summary]

    Args:
        r0 (float): horizontal distance from center of hemisphere
        R (float): radius of lens
        n1 (float): index of refraction of lens
    """
    
    # calculate incidence angle on lenslet 
    # we are assuming rays come vertically down
    theta0 = asin(r0/R)
    
    # Snell's law
    theta1 = asin( sin(theta0)/ n1)
    
    # calculate beta
    beta = pi/2 - theta1 + theta0
    
    # calculate location of ray from the center of hemisphere
    r = R * sin(theta1) / sin(beta)
    
    # outputs
    #print("Theta 0: %.1f " %(degrees(theta0)) )
    #print("Theta 1: %.1f " %(degrees(theta1)) )
    #print("r0: %4.2f  r: %4.2f " %(r0, r))
    
    return r


def ray_path_single_ar(r0:float, R:float, t:float, n1:float, n2:float):
    """This code approximates the path of a ray thorugh a lenslet coated with an AR

    Args:
        r0 (float): horizontal distance of incoming ray from center of hemisphere
        R (float): radius of lens
        t (float): thickness of AR
        n1 (float): index of refraction of lens
        n2 (float): index of refraction of AR
    """
    
    # calculate incident angle on AR
    theta0 = asin(r0/(R+t))
    
    # snells law at the vacuum/ar interface
    theta02 = asin(sin(theta0)/n2)
    
    # calculate *approximate* distance covered by ray through the AR
    l = t/cos(theta02)
    
    # *approximate* theta21 as theta02, where theta21 is the incident angle on the lenslet.
    theta21 = theta02
    
    # snell's law at the lenslet/ar interface
    theta1 = asin(n2/n1 * sin(theta21))
    
    # Now we need to find the distance from incident angle to center of hemisphere again (r1)
    phi = theta0 - theta02
    y = l * sin(phi)
    
    r1 = r0 - y
    
    # finally, we can find r
    psi = asin(r1/R)
    beta = pi/2 - theta1 + psi
    
    r = R * sin(theta1)/sin(beta)
    
    return r
    
    
    
    

if __name__ == "__main__":
    main()