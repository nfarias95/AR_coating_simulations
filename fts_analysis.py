# code to get transmission data from fts
import numpy as np
from matplotlib import pyplot as plt
import csv


DATA_FOLDER = 'C:/Users/nicol/Documents/00Research/Data/Thermal_Sprayed_Samples/2022_march_31_two_layer_coupons/'

LOAD_PREDICTION = True
pred_file_name = 'FTSprediction.txt'

def main():
    print("\nanalyzing fts data...")
    
    cal_files = np.array(['1calib.fft', '2calib.fft', '5calib.fft', '9calib.fft']) # calibration (open, no samples)
    sample_files = np.array(['3coupon.fft', '4coupon.fft', '7coupon.fft', '8coupon.fft']) # sample
    
    num_samples = len(sample_files)
    print("number of samples: ", num_samples)
    
    # load the data
    trans_sum = 0 # to take the average
    for i in range(0, num_samples):
        cal_freqs, cal_amp = load_data(cal_files[i])
        sample_freqs, sample_amp = load_data(sample_files[i])
        
        # calculate the transmission
        trans = sample_amp/cal_amp
        trans_sum = trans_sum + trans
    
    trans_mean = trans_sum / num_samples
    
    plt.scatter(cal_freqs/1e9, trans_mean, label='FTS data', s=8)
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Transmission")
    plt.xlim([30, 180])
    plt.ylim([0, 1.2])
    plt.grid()
    
    if LOAD_PREDICTION:
        pred_freqs, pred_amp = load_prediction_data(pred_file_name)
        plt.plot(pred_freqs, pred_amp, label='prediction', color='red')
        plt.legend()

    
    
    plt.show()
    print("the end. \n\n")
    
def load_prediction_data(pred_file_name : str):
    """_summary_

    Args:
        pred_file_name (str): _description_
    """
    freqs = []
    amplitudes = []
    
    with open(DATA_FOLDER + pred_file_name) as file:
        csvreader = csv.reader(file, delimiter=" ")
        i=-1
        for line in csvreader:
            i+=1
            if i >1:
                freqs.append(float(line[-3]))
                amplitudes.append(float(line[-1]))
    
    return np.array(freqs), np.array(amplitudes)
    

def load_data(filename):
    # load the .if data of a run from comma separated value text file
    file = open(DATA_FOLDER + filename)     
    
    csvreader = csv.reader(file)
    
    freqs = []
    amplitudes = []
    
    for row in csvreader:
        freqs.append(float(row[0]))
        complex_value = complex(row[2])
        amplitudes.append(np.absolute(complex_value))

    file.close()
    
    return np.array(freqs), np.array(amplitudes)

if __name__ == "__main__":
    main()