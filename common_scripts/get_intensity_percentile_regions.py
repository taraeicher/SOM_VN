import region_defs
import pickle as pkl
import numpy as np
import math
from tqdm import tqdm
import sys

def main():
    # Input
    region_file = pkl.load(open(sys.argv[1], "rb"))
    percentile = float(sys.argv[2])
    percentile_cutoff_file = open(sys.argv[3], "w")
    
    # Percentile
    all_signal = concatenate(region_file)
    percentile_cutoff = get_intensity_percentile(percentile, all_signal)
    percentile_cutoff_file.write(str(percentile_cutoff))
    percentile_cutoff_file.close()

"""
Concatenate all signal by accessing the "signals" field in the
region definitions.
"""
def concatenate(region_list):
    all_signal = []
    print("Extending list of all signal intensities")
    for region in tqdm(region_list):
        all_signal.extend(list(region.signals))
    
    return all_signal

"""
Returns the Xth percentile of intensity for all records in the file. 
NOTE: 1000000 is the highest possible RPKM intensity.
"""
def get_intensity_percentile(percentile, signals):
    fine_bin_count = 4
    max_threshold = 1000000 * fine_bin_count
    counts = np.zeros(max_threshold)        
    file_line_count = 0
    
    print("Obtaining signal percentile")
    for signal in tqdm(signals):
        
            #Increment the count of bins.
            file_line_count += 1
            
            #Increment the appropriate location by the number of bins.
            bin = int(math.ceil(signal))
            counts[bin] += 1

    #Find percentile of maxes.
    target_count = np.sum(counts[1:]) * percentile
    running_sum = 0
    i = 1
    percentile_found = False
    max_sig_percentile = 0
    
    while i < len(counts) - 1 and not percentile_found:
        running_sum += counts[i]
        if running_sum >= target_count:
            max_sig_percentile = i
            percentile_found = True
        i += 1
        
    #Rewind file and return values.
    retval = max_sig_percentile / fine_bin_count
    return retval
    
if __name__ == "__main__":
    main()