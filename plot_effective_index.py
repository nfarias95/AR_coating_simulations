# Code to plot the effective index of refraction
import numpy as np
import math
from matplotlib import pyplot as plt

def main():
    FONT_SIZE = 20
    plt.rcParams.update({"font.family": "Times New Roman"})
    plt.rcParams.update({'font.size': FONT_SIZE})
    fig = plt.figure(figsize=(8,6))


    um = "\u03BCm"
    
    print("Let's calculate the effective index of refraction of a metamaterial given a shape")
    
    r_array = np.linspace(30, 80, 6)
    print("R: ", r_array)
    
    #neff_270 = np.array([1.815 ,1.662, 1.61 , 1.508, 1.508, 1.456])
    #neff_475 = np.array([1.86, 1.764, 1.713, 1.662, 1.661,  1.61 ])
    #neff_680 = np.array([1.918, 1.815, 1.713, 1.662, 1.508, 1.456 ])
    
    neff_270 = [1.862, 1.758, 1.656, 1.600, 1.559, 1.532]
    neff_475 = [1.869, 1.765, 1.676, 1.621, 1.586, 1.545]
    neff_680 = [1.840, 1.758, 1.682, 1.621, 1.593, 1.552]
    
    plt.scatter(r_array, neff_270, marker="*", color="k", s=100, label="270 "+um)
    plt.scatter(r_array, neff_475, marker="^", s=50, label="475 "+um)
    plt.scatter(r_array, neff_680, label="680 "+um)
    plt.legend()
    plt.xlabel("Cylinder radius [\u03BCm] ")
    plt.ylabel("Effective index of refraction")
    plt.show()
    
    print("the end.")
    
if __name__ == "__main__":
    main()