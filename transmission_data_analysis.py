# This code analyzes transmission data from HFSS
import csv
import numpy as np
from matplotlib import pyplot as plt


# MAIN
def main():
    
    print("Welcome. We'll read transmission data")
    
    # insert file path here:
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/July2021/"
    # insert file name here:
    f_name = "date_07052021_corrected_1cyl_s_10_tr_var_th_680.csv"
    # complete file location:
    f_loc = f_path + f_name
    
    #insert desired column here:
    freq_column = 2
    data_column = 3
    
    #insert desired frequency band here
    #band = np.array([34, 99])  # FB1 
    #band = np.array([60, 162]) # FB2
    band = np.array([77, 224]) # FB3
    
    #read the data
    freq_array, trans_array, title = read_csv_file(f_loc, freq_column, data_column)
       
    print("Frequency: \n", freq_array)
    print("Transmission: \n", trans_array)
       
    # Plot the data
    plot_results(freq_array, trans_array, title, band)
    
    # Calculate statistics
    mean_transmission = np.mean( trans_array[band[0]: band[1] ])
    print("Mean transmission: %.1f " % (mean_transmission * 100) )
    
    print("The end.")
    
    
########### end of main
    
def plot_results(freq_array, trans_array, title, band):
    
    plt.scatter(freq_array, trans_array)
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Transmission")
    plt.title(title)
    
    # add a line at 95%
    plt.axhline(y=0.95, color = 'r', linestyle='-.', linewidth = 0.5)
    
    # add lines at the band edges
    plt.axvline(x=band[0], color='k', linestyle='-.', linewidth = 0.5)
    plt.axvline(x=band[1], color='k', linestyle='-.', linewidth = 0.5)
        
    plt.show()
    
def read_csv_file(f_loc:str, x_column:int, y_column:int):
    """ This function reads a csv file

    Args:
        f_loc ([type]): file location
        freq_column (int): column in data file where independent variable is located
        data_column (int): column where dependent variable is located

    Returns:
        [type]: [description]
    """
    
    
    fields = []
    rows = []
    
    with open(f_loc, 'r') as csvfile:
        print("Reading data")
        #creating a csv reader object
        csvreader = csv.reader(csvfile)
        
        #extracting field names through first row
        fields = next(csvreader)
        
        #extracting each data row one by one
        for row in csvreader:
            rows.append(row)
            
        #get total number of rows
        print("\nTotal number of rows: %d " % (csvreader.line_num))
        
    # printing th field names
    print('\nField names are: ' + ', ' .join(field for field in fields))
    
    # printing first 5 rows
    print("\n First 5 rows are: \n")
    counter = -1
    for row in rows[:5]:
       counter = counter + 1
       print(rows[counter])
       #for col in row:
           #print("%10s " % col)
       print('\n')
        
        
    freq_array = np.zeros(csvreader.line_num -1) # frequency array
    trans_array = np.zeros(csvreader.line_num -1) # transmission array
    
    # put frequency and transmission into separate arrays
    title = fields[y_column]
    counter = -1
    for row in rows:
        counter = counter + 1
        freq_array[counter] = float(row[x_column])
        trans_array[counter] = float(row[y_column])
    
    
    
    
    return freq_array, trans_array, title
        
    
    
if __name__ == "__main__":
    main()
    
    