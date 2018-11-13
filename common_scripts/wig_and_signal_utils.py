import numpy as np
import math
import sys

"""
Estimate the number of peaks in the region.
Do this by checking the number of trough-peak-trough triples,
where trough and peak are signals within a threshold of the
max and min.
""" 
def find_peak_count(region, threshold, lowest):
    #cutoffs = [0.5, 1.25]
    #max_val = np.max(np.asarray(region))
    #min_val = np.min(np.asarray(region))
                
    at_trough = True
    peak_count = 0
    i = 0
    
    #If there is not a significant difference between max and min, we don't worry about peaks.
    #if max_val >= threshold * min_val and max_val - threshold >= min_val:
    for sig in region:
        #If the last thing we saw was a trough and we found a peak, update the last thing
        #we saw to a peak.
        if at_trough and sig >= threshold:#cutoffs[0] * threshold:
            at_trough = False
            
        #If the last thing we saw was a peak and we found a trough, update the last thing
        #we saw to a trough and update the peak count.
        if (not at_trough) and sig < threshold:#cutoffs[1] * min_val:
            at_trough = True
            peak_count += 1
        i += 1
        #If we did not find a trough-peak-trough combo, set number of peaks to the fraction of the threshold.
        if peak_count == 0:
            peak_count = (np.max(region) - lowest) / (threshold - lowest)
        return peak_count
    
# """
# Returns the Xth percentile of intensity for all records in the file. 
 # """
# def get_intensity_percentile(percentile, file):
    
    # #Count of the number of inputs with each FDI and max.
    # #Note: for window size of 10, fdi_threshold is about 6.
    # #FDI scaling factor allows us to store the FDI data at
    # #a finer granularity.
    # max_threshold = 10000
    # counts = np.zeros(max_threshold * 4)        
    # file_line_count = 0
    
    # #Use each entry in the file to calculate running metadata.
    # next_line = file.readline()
    
    # while next_line:
        # split_line = next_line.split()
        # if len(split_line) == 2:
            # val = float(split_line[1])
            
            # #Increment the count of lines in the file.
            # file_line_count += 1
            
            # #Get the maximum value and increment the appropriate location in the array.
            # counts[int(math.ceil(4 * val))] += 1
            
            # #Let the user know we are still processing.
            # if file_line_count % 10000 == 0:
                # sys.stdout.write('.')
                
        # #Read the next line.
        # next_line = file.readline()
        
    # #Find percentile of maxes.
    # target_count = int(file_line_count * percentile)
    # running_sum = 0
    # i = 0
    # percentile_found = False
    # max_sig_percentile = 0.0
    # while i < len(counts) - 1 and not percentile_found:
        # running_sum += counts[i]
        # if running_sum >= target_count:
            # max_sig_percentile = i / 4.0
            # percentile_found = True
        # i += 1
    # i = 0
    # percentile_found = False
    # running_sum = 0
        
    # #Return values.
    # #print(max_sig_percentile)
    # return max_sig_percentile
    
"""
Returns the Xth percentile of intensity for all records in the file. 
 """
def get_intensity_percentile(percentile, file, lowest, fine=False):
    
    #Count of the number of inputs with each FDI and max.
    #Note: for window size of 10, fdi_threshold is about 6.
    #FDI scaling factor allows us to store the FDI data at
    #a finer granularity.
    fine_bin_count = 4
    max_threshold = 10000
    to_subtract = lowest
    if fine:
        max_threshold *= fine_bin_count
        to_subtract *= fine_bin_count
    counts = np.zeros(max_threshold)        
    file_line_count = 0
    
    #Use each entry in the file to calculate running metadata.
    next_line = file.readline()
    
    while next_line:
        split_line = next_line.split()
        if len(split_line) == 2:
            val = float(split_line[1])
            
            #Increment the count of lines in the file.
            file_line_count += 1
            
            #Get the maximum value and increment the appropriate location in the array.
            if fine:
                val *= fine_bin_count
            bin = int(math.ceil(val - to_subtract))
            counts[bin] += 1
            
            #Let the user know we are still processing.
            if file_line_count % 10000 == 0:
                sys.stdout.write('.')
                
        #Read the next line.
        next_line = file.readline()
        
    #Find percentile of maxes.
    target_count = int(file_line_count * percentile)
    running_sum = 0
    i = 0
    percentile_found = False
    max_sig_percentile = 0
    
    while i < len(counts) - 1 and not percentile_found:
        running_sum += counts[i]
        if running_sum >= target_count:
            max_sig_percentile = i
            percentile_found = True
        i += 1
    i = 0
    percentile_found = False
    running_sum = 0
        
    #Rewind file and return values.
    file.seek(0)
    retval = max_sig_percentile + to_subtract
    if fine:
        retval /= 4
    return retval
    
#Get the cross-correlation metric between the two clusters.
def get_crosscorr(clust1, clust2, delay, threshold, max_cutoff, use_max_cutoff, two_way, minimum):
    
    #Get the section of the cluster that does not include the delay.
    clust1_array = np.asarray(clust1)
    clust2_array = np.asarray(clust2)

    #If the clusters are the same length, use subarrays to simulate the shift.
    #Otherwise, move the smaller array across the larger one.
    if two_way:
        if delay < 0:
            #Simulate the shift by moving the second cluster to the right.
            clust1_array = np.asarray(clust1[0:(len(clust1) + delay)])
            clust2_array = np.asarray(clust2[abs(delay):len(clust2)])
            
        else:
            #Simulate the shift by moving the first cluster to the right.
            clust2_array = np.asarray(clust2[0:(len(clust2) - delay)])
            clust1_array = np.asarray(clust1[delay:len(clust1)])
    elif len(clust1_array) > len(clust2_array):
        if delay < 0:
            raise ValueError("Cannot use a negative delay")
            
        else:
            clust1_array = np.asarray(clust1[delay:(len(clust2) + delay)])
    else:
        if delay < 0:
            raise ValueError("Cannot use a negative delay")
            
        else:
            clust2_array = np.asarray(clust2[delay:(len(clust1) + delay)])
    
    #Only calculate the cross-correlation if the sub-regions of both clusters contain the max
    #and the maximums are within the threshold.
    #Else, throw an error.
    both_contain_max = np.max(np.asarray(clust1)) == np.max(clust1_array) and np.max(np.asarray(clust2)) == np.max(clust2_array)
    max_of_two = max(np.max(np.asarray(clust1)), np.max(np.asarray(clust2)))
    min_of_two = min(np.max(np.asarray(clust1)), np.max(np.asarray(clust2)))
    maxes_within_threshold =  min_of_two / max_of_two > threshold or max_of_two == minimum

    if max_of_two < max_cutoff and use_max_cutoff:
        return 1
    elif both_contain_max and maxes_within_threshold:
        #Calculate the cross-correlation pieces.    
        numerator = np.sum(clust1_array * clust2_array) - np.sum(clust1_array) * np.sum(clust2_array) / len(clust2_array)
        denominator_1 = np.sum(np.square(clust1_array)) - np.square(np.sum(clust1_array)) / len(clust1_array)
        denominator_2 = np.sum(np.square(clust2_array)) - np.square(np.sum(clust2_array)) / len(clust2_array)
        denominator = 1
        if denominator_1 > 0 and denominator_2 > 0:
            denominator = np.sqrt(denominator_1 * denominator_2)
        else:
            numerator = 0
    else:
        raise ValueError("Maximum signal intensity not contained in region")
        
    #Return the cross-correlation result.
    return numerator / denominator