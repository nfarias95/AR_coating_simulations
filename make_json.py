import json

outputfile_name = "date03142022_single_layer_cylinder_td_160_h_680"
outputfile_path = "C:/Users/nicol/Documents/00Research/Data/MetamaterialSims/March2022 Effective Index/"
outputfile_type = ".txt"
outputfile_loc = outputfile_path + outputfile_name + outputfile_type


x = {
    "filename" : outputfile_name,
    "path" : outputfile_path, 
    "data_label": "TD=160, H=680",
    "freq_column" : 2,
    "data_column" : 3,
    "var_column" : -1,
    "var_label" : "none",
    "f_convert" : 1e9,  # given frequency in the data is in GHz, but we'll use Hz for our calculations
    "freq_min": -1,  # will be the minimum frequency of the data given (no need to trim)
    "freq_max": -1, # will be the maximum of the data given (no need to trim)
    # SETUP PARAMETERS
    'inc_angle': 0, # degrees
    # NOW INSERT THE CHARACTERISTICS OF THE MATERIAL
    'opt_type' : 'lenslet',
    'n_substrate' : 3.4, # index of refraction
    'lt_substrate' : 0.00005,
    't_substrate' : 9999,  # [meters] (doesn't matter for lenslet)
    'n_ar_1' : 1.7, #  index of refraction of coating (roughly)
    'n_ar_2' : -1, # not applicable
    't_ar_1' : 680e-6, # thickness of coating #1
    't_ar_2' : -1 # n/a
    
    
}


with open(outputfile_loc, 'w') as outfile:
    json.dump(x, outfile)

print("The end. \n\n")