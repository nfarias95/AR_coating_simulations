# This code simulates the performance of anti-reflection coating

import numpy as np
import matplotlib.pyplot as plt
import math


# source: https://farside.ph.utexas.edu/teaching/em/lectures/node104.html#:~:text=is%20a%20unit%20vector%20pointing%20in%20the%20direction%20of%20wave%20propagation.&text=Let%20us%20investigate%20what%20happens,this%20boundary%20from%20medium%201.

def main():
    
    print("Hello, let's calculate transmissions at an interface between vacuum and a dielectric")
    
    # index of refraction
    n1 = 1
    n2 = 3.4 # silicon
    
    # incidence angle
    thetas_i = np.array([0, 15, 30, 45, 60, 75, 90]) * np.pi/180


    # calculate transmission at each angle
    print("\n\n Theoretical Transmission (Maxwell Eqns)")
    print("          Angle       Transmission")
    for i in range(len(thetas_i)):
        theta_i = thetas_i[i]
        # snells law: sin(theta1) * n1 = sin(theta2) * n2
        theta_t = np.arcsin( n1/n2 * np.sin(theta_i))
        # calculate alpha and beta
        beta = n2/n1
        alpha = np.cos(theta_t) / np.cos(theta_i)
        
        T = alpha * beta * (2/(1+alpha*beta))**2 # equation 1262
    
        print("            %4.2f   %4.2f " %(theta_i * 180/np.pi, T) )
    
    print("The end. \n\n")

if __name__ == "__main__":
    main()
    





