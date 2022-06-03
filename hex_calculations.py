import math
import numpy as np

pi = np.pi
sqrt = np.sqrt

def main():
    print("\nhi")
    
    r = 30 # radius of circle
    s = 10 # space between holes
    
    a_hole = 25 # length of side of hexagonal hole
    
    fc = f_circle(r, s)
    
    fhex = f_hexagon(a_hole, s)
    
    print("f circle: ", fc)
    print("f hex: ", fhex)
    

def f_hexagon(a_hole, s):
    
    Ahole = 3*sqrt(3)/2 * a_hole**2
    
    h_hole = sqrt(3)/2 * a_hole
    
    h_cell = h_hole + s/2
    
    a_cell = 2/sqrt(3) * h_cell
    
    Acell = 3*sqrt(3)/2 * a_cell**2
    
    return Ahole/Acell

def f_circle(r, s):
    
    Acircle = pi * r**2
    
    h = r + s/2
    a = 2 / sqrt(3) * h
    
    Ahex = 3 * sqrt(3)/2 * a **2
    
    f = Acircle/Ahex
    return f

if __name__ == "__main__":
    main()
    