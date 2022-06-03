# This code analyzes transmission data from HFSS
# It allows you to compare transmissions with different parameters (e.g., thickness of an AR)
import csv
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from useful_functions import *

# MAIN
SAVE_DATA = False

# PLOT PARAMETERS
FONT_SIZE = 20
LEGEND_SIZE = 18
PLOT_95_LINE = True

# set plot font
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': FONT_SIZE})
plt.rc('legend', fontsize=LEGEND_SIZE)    # legend fontsize
#fig = plt.figure(figsize=(8,6))
um = "\u03BCm"

def main():
    print("Welcome. We'll read transmission data")
    
    # FILE LOCATION AND DESCRIPTION    
    f_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/October_2021_redoing_conical_sims/"
    f_name = "date03042022_conical_TR_var_TH_1000_s_10_selected_variations.csv"
    labels = ["80 " + um, "120 " + um, "160 " + um, "280 " + um, "320 " + um, "400 " + um]
    #plot_title = "Effect of changing cone radius"
    
    #Columns from csv file:
    freq_column = 4 # column containing frequency data
    data_column = 5 # column containing the desired data (transmission) data
    var_column = 3 # column containing a variation parameter (ex: radius, height, spacing). For simulations only. Set -1 if no parameters are varying
    
    # GET COMPLETE FILE LOCATION:
    f_loc = f_path + f_name
    
    #insert desired frequency band here
    #band = np.array([34, 99])  # FB1 
    #band = np.array([60, 162]) # FB2
    #band = np.array([77, 224]) # FB3
    band = np.array([50 , 250])
    
    #read the data
    freq_array, trans_array, var_array, title = read_csv_file(f_loc, freq_column,  data_column, var_column)
       
    # split the data if there are varying parameters
    
    if var_column >= 0:
         freq_array, trans_array, var_values  = split_data(var_array, freq_array, trans_array)
    else:
        var_values = []
         
    
    
    #Simple plot
    plt.figure(9, figsize=(8,6.5))
    #marker_type = ["."]
    marker_type = ["x", ".", "o", "d", "v", "*", ">"]
    for i in range(0, len(var_values)):
            plt.plot(freq_array, trans_array[:, i], label=labels[i], marker=marker_type[i], markersize=5)
            # save the data if you want
            if SAVE_DATA:
                print("fpath + labels[i]", f_path + labels[i])
                output_file_name = f_name[0:-4] + "_" + labels[i] + ".csv"
                write_data_to_csv_files(freq_array, trans_array[:,i], "freq [ghz]", "transmission",  output_file_name, f_path )
            
    # Plot line at 95%
    if PLOT_95_LINE:
        xmin = 50
        xmax = 250
        x95 = np.linspace( xmin, xmax, 3)
        y95 = 0.95 * np.ones(3)
        #plt.plot(x95, y95, '--', markersize=0.00001, color='k')  
        #plt.axhline(y=0.95, color = 'k', linestyle='--', linewidth = 0.5)
        plt.axhline(y=0.95, color = 'k', linestyle='--', markersize=0.00001)
        
        plt.xlim([xmin, xmax])
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Transmittance")
    #plt.title(plot_title,  size='15')
    plt.legend(loc = 'lower right')
    plt.show()
    
    # Calculate statistics
    mean_transmission = np.mean( trans_array[band[0]: band[1] ])
    print("Mean transmission: %.1f " % (mean_transmission * 100) )
    
    print("The end.")
    
    
########### end of main
def split_data(var_array:np.ndarray, freq_array:np.ndarray, trans_array:np.ndarray):
    """[summary]

    Args:
        var_array (np.ndarray): array with varying parameter
        freq_array (np.ndarray): array with all frequency data
        trans_array (np.ndarray): array with all transmission data
    """
    print("Separating the data.")
    numpoints = len(var_array) # total number of points in var_array
    
    var_old = var_array[0]
    
    var_values = np.array([var_old]) # start an array to store the varying values
    
    for i in range(0, numpoints):
        var = var_array[i]
        
        if var != var_old:
            #print("var: ", var)
            var_values = np.append(var_values, var)
            
        
        var_old = var
    print("Varying values: ", var_values)
    #calculate number of variations:
    num_var = len(var_values)
    num_new_rows = int( len(var_array)/num_var)
    print("Number of varitations: ", num_var)
    print("Number of rows: ", num_new_rows)
    
    # divide the data into more arrays
    freq_array_new = freq_array[0:num_new_rows] # the frequency array does not change, take the first ones
    trans_array_new = np.zeros((num_new_rows, num_var))
    
    # send values over
    for i in range(0, num_var):
        trans_array_new[:, i] = trans_array[i*num_new_rows:(i+1)*num_new_rows]
        
    return freq_array_new, trans_array_new, var_values    
    
    
    
def plot_results(freq_array, trans_array, var_values, title, band):
    
    plt.figure(10, figsize=(8,6))
    if len(var_values) > 0:
        for i in range(0, len(var_values)):
            plt.scatter(freq_array, trans_array[:, i], label=str(var_values[i]))
    else:
        plt.scatter(freq_array, trans_array)
        
    plt.xlabel("Frequency [GHz]", fontname="Times New Roman")
    plt.ylabel("Transmission", fontname="Times New Roman")
    #plt.title(title, fontname="Times New Roman")
    plt.legend()
    
    # add a line at 95%
    plt.axhline(y=0.95, color = 'r', linestyle='-.', linewidth = 0.5)
    
    # add lines at the band edges
    plt.axvline(x=band[0], color='k', linestyle='-.', linewidth = 0.5)
    plt.axvline(x=band[1], color='k', linestyle='-.', linewidth = 0.5)
        
    plt.show()
    
    
def read_csv_file(f_loc:str, x_column:int, y_column:int, z_column:int):
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
        #print("Reading data")
        #creating a csv reader object
        csvreader = csv.reader(csvfile)
        
        #extracting field names through first row
        fields = next(csvreader)
        
        #extracting each data row one by one
        for row in csvreader:
            rows.append(row)
            
        #get total number of rows
        #print("\nTotal number of rows: %d " % (csvreader.line_num))
        
    # printing th field names
    #print('\nField names are: ' + ', ' .join(field for field in fields))
    
    # printing first 5 rows
    #print("\n First 5 rows are: \n")
    # counter = -1
    # for row in rows[:5]:
    #    counter = counter + 1
    #    print(rows[counter])
    #    #for col in row:
    #        #print("%10s " % col)
    #    print('\n')
        
        
    freq_array = np.zeros(csvreader.line_num -1) # frequency array
    trans_array = np.zeros(csvreader.line_num -1) # transmission array
    var_array = np.zeros(csvreader.line_num -1) # varying field (ex., radius, height, spacing)
    
    # put frequency and transmission into separate arrays
    title = fields[y_column]
    counter = -1
    for row in rows:
        counter = counter + 1
        freq_array[counter] = float(row[x_column])
        trans_array[counter] = float(row[y_column])
        if z_column >= 0:
            var_array[counter] = float(row[z_column])
    
    
    
    return freq_array, trans_array, var_array, title
        
    
    
if __name__ == "__main__":
    main()
    
    