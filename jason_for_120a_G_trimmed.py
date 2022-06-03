# This program creates Jason file for the data we want
import json

outputfile_name = "120a_G_trimmed"
outputfile_path = "C:/Users/nicol/Documents/00Research/Data/Thermal_Sprayed_Samples/2021_august_CTS_Hayward_first_batch/"
outputfile_type = ".txt"
outputfile_loc = outputfile_path + outputfile_name + outputfile_type

# create a dictionary
x = {
    "filename" : outputfile_name,
    "path" : outputfile_path, 
    "data_label": "120a",
    "freq_column" : 0,
    "data_column" : 4,
    "var_column" : -1,
    "var_label" : "none",
    "f_convert" : 1e9,  # given frequency in the data is in GHz, but we'll use Hz for our calculations
    "freq_min": -1,  # will be the minimum frequency of the data given (no need to trim)
    "freq_max": -1, # will be the maximum of the data given (no need to trim)
    # SETUP PARAMETERS
    'inc_angle': 5.16, # degrees
    # NOW INSERT THE CHARACTERISTICS OF THE MATERIAL
    'opt_type' : 'coupon',
    'n_substrate' : 3.11, # index of refraction of alumina
    'lt_substrate' : 0.02,
    't_substrate' : 0.00635,  # [meters]
    'n_ar_1' : 2.378, #  index of refraction of coating (roughly)
    'n_ar_2' : -1, # not applicable
    't_ar_1' : 0.000490728, # thickness of coating #1
    't_ar_2' : -1 # n/a
    
    
}

#y = json.dumps(x)

with open(outputfile_loc, 'w') as outfile:
    json.dump(x, outfile)

print("The end. \n\n")