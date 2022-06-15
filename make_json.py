import json

outputfile_name = "example"
outputfile_path = "C:/Users/nicol/Documents/00Research/PythonCode/AR_coating_simulations/"
outputfile_type = ".txt"
outputfile_loc = outputfile_path + outputfile_name + outputfile_type


x = {
    "filename" : outputfile_name, # should be the same as the csv data file name
    "path" : outputfile_path,  # should be the same as the path to the csv data 
    "data_label": "TD=160, H=475", # a description to your data
    "freq_column" : 2, # column where frequency data is located [count from 0]
    "data_column" : 3, # column where data (transmission) is located [count from 0]
    "var_column" : -1, # column where a varying parameter (for example, the thickness of a layer) is changing. if none, use -1
    "var_label" : "none", # label of the parameter that is varying. write "none" if not applicable
    "f_convert" : 1e9,  # given frequency in the data is in GHz, but we'll use Hz for our calculations
    "freq_min": -1,  # will be the minimum frequency of the data given (-1 if no need to trim)
    "freq_max": -1, # will be the maximum of the data given (-1 if no need to trim)
    # SETUP PARAMETERS
    'inc_angle': 0, # degrees
    # NOW INSERT THE CHARACTERISTICS OF THE MATERIAL
    'opt_type' : 'lenslet', # can be lenslet, lens or (???)
    'n_substrate' : 3.4, # index of refraction of substrate. 3.4 for silicon, 3.1 for alumina
    'lt_substrate' : 0.00005, # loss tangent of substrate.
    't_substrate' : 9999,  # thickness of substrate [meters] (doesn't matter for lenslet, since the silicon is assumed to be infinitely thick in one direction)
    'n_ar_1' : 1.7, #  approximate index of refraction of coating (actual value can be found with fit)
    'n_ar_2' : -1, # approximate index of refraction of second coating (if single layer, write -1)
    't_ar_1' : 475e-6, # thickness of coating [m] #1 MAKE THIS THE SAME AS YOUR METAMATERIAL HEIGHT (e.g., hole depth)
    't_ar_2' : -1, # thickness of AR coating [m] #2 (n/a if single layer, write -1)
    # FIT PARAMETERS
    'tol':  0.2, # tolerance: how far we can be from original guess. range: 0-1
    'N':  200 # number of data points to be analyzed per fit. the more, the more accurate and also the slower!
    
}


with open(outputfile_loc, 'w') as outfile:
    json.dump(x, outfile)

print("Created json file: ", outputfile_loc)

print("The end. \n\n")