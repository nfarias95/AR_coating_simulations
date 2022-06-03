# file containing useful functions

import csv
import numpy as np
from matplotlib import pyplot as plt

def write_data_to_csv_files(x_array, y_array, xname, yname, fname, folder):
    fileloc=folder + fname
    print("Saving data...")
    print("destination file name: ", fileloc)
    
    if len(x_array) != len(y_array):
        print("\n\nX and Y should have the same number of elements!\n\n")
    
    npoints = len(x_array)
    
    
    with open(fileloc, 'w', newline='') as f:
        #create the csv writer
        writer = csv.writer(f)
        
        # write the rows
        row = [xname, yname]
        writer.writerow(row)
        for i in range(0, npoints):
            row = [x_array[i], y_array[i]]
            #print(row)
            writer.writerow(row)
            
def plot_95_per_cent_line():
    xmin = 0
    xmax = 250
    x95 = np.linspace( xmin, xmax, 3)
    y95 = 0.95 * np.ones(3)
    plt.plot(x95, y95, '--', markersize=0.00001, color='k') 
    plt.xlim([xmin, xmax])
    
    
def shade_bands(ax, f1:float, f2:float, clr:str):
    """ Function to add a shade between frequency bands in a plot

    Args:
        i (float): a figure counter
        f1 (float): frequency start [GHz]
        f2 (float): frequency band [GHz]
        clr (str): color of shaded region
    """
    
    X = np.linspace(f1, f2, 3)
    Y = np.array([1, 1, 1])
    ax.fill_between(X, 0, Y, color=clr, alpha=0.2)