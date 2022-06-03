import numpy as np
from matplotlib import pyplot as plt

def my_main():
    
    stan_dist = [1.5, 2, 2.5, 3, 3.5, 4] # x
    
    feed_rate = [2, 
                 3, 
                 4, 
                 5, 
                 6 ] # y
    
    eps = [ 
            [0,   2.79, 2.93, 3.06, 3.11, 3.11], #2
            [2.6, 2.51, 2.90, 0,    2.84, 2.87], #3
            [0,   2.78, 2.93, 0,    0,    0] , #4
            [0,   2.7,  2.9,  3.07, 0,    0],#5
            [0,   2.39, 2.36, 2.36, 0,    0]
        ]
    
    
    plt.pcolormesh(eps, cmap="magma")
    plt.colorbar()
    
    x_tick_location = np.linspace(0.5, 5.5, 6)
    y_tick_location = np.linspace(0.5, 4.5, 5)
    
    plt.xticks(x_tick_location, stan_dist)
    plt.yticks(y_tick_location, feed_rate)
    
    plt.xlabel("Standoff Distance [in]")
    plt.ylabel("Powder Feed Rate [lb/hr]")
    plt.title("Fillite Dielectric Constant")
    
    # ----------
    plt.show()
    print("the end")
    
if __name__ == "__main__":
    my_main()