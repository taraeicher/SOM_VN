import numpy as np
import math
import sys

"""
Estimate the number of crests in the region.
Do this by checking the number of trough-crest-trough triples,
where trough and crest are signals within a threshold of the
max and min.
""" 
def find_crossing_count(region, threshold):
               
    is_below = True
    crossing_count = 0

    for sig in np.nditer(region):
        #If the last thing we saw was a trough and we found a crest, update the last thing
        #we saw to a crest.
        if is_below and sig >= threshold:
            is_below = False

        #If the last thing we saw was a crest and we found a trough, update the last thing
        #we saw to a trough and update the crest count.
        if (not is_below) and sig < threshold:
            is_below = True
            crossing_count += 1

    #If we did not find a trough-crest-trough combo, set number of crests to the fraction of the threshold.
    if crossing_count == 0:
        crossing_count = min(1, np.max(region) / threshold)
        
    return crossing_count
        
"""
Returns the Xth percentile of intensity for all records in the file. 
NOTE: 1000000 is the highest possible RPKM intensity.
"""
def get_intensity_percentile(percentile, file):
   
    fine_bin_count = 4
    max_threshold = 1000000 * fine_bin_count
    counts = np.zeros(max_threshold)        
    file_line_count = 0
    
    #Use each entry in the file to calculate running metadata.
    next_line = file.readline()
    
    while next_line:
        split_line = next_line.split()
        if len(split_line) == 4:
            val = float(split_line[3]) * fine_bin_count

            #Increment the count of lines in the file.
            file_line_count += 1
            
            #Get the maximum value and increment the appropriate location in the array.
            bin = int(math.ceil(val))
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
    retval = max_sig_percentile / fine_bin_count
    return retval
    
#Get the cross-correlation metric between the two clusters.
def get_crosscorr(shape1, shape2, delay, threshold, max_ratio_cutoff, two_way, percentile_cutoff):

    #If the clusters are the same length, use subarrays to simulate the shift.
    #Otherwise, move the smaller array across the larger one.
    shape1_subset = shape1[0:len(shape1)]
    shape2_subset = shape2[0:len(shape2)]
    if two_way:
        if delay < 0:
            #Simulate the shift by moving the second cluster to the right.
            shape1_subset = shape1[0:(len(shape1) + delay)]
            shape2_subset = shape2[abs(delay):len(shape2)]
            
        else:
            #Simulate the shift by moving the first cluster to the right.
            shape2_subset = shape2[0:(len(shape2) - delay)]
            shape1_subset = shape1[delay:len(shape1)]
    elif len(shape1_subset) > len(shape2_subset):
        if delay < 0:
            raise ValueError("Cannot use a negative delay")
            
        else:
            shape1_subset = np.asarray(shape1[delay:(len(shape2) + delay)])
    else:
        if delay < 0:
            raise ValueError("Cannot use a negative delay")
            
        else:
            shape2_subset = np.asarray(shape2[delay:(len(shape1) + delay)])
    
    #Only calculate the cross-correlation if the sub-regions of both clusters contain the max
    #and the maximums are within the threshold.
    #Else, throw an error.
    both_contain_max = np.max(shape1) == np.max(shape1_subset) and np.max(shape2) == np.max(shape2_subset)
    max_of_two = max(np.max(shape1), np.max(shape2))
    min_of_two = min(np.max(shape1), np.max(shape2))
    maxes_within_threshold =  min_of_two / max_of_two > max_ratio_cutoff or max_of_two == percentile_cutoff

    elif both_contain_max and maxes_within_threshold:
        #Calculate the cross-correlation pieces.    
        numerator = np.sum(shape1_subset * shape2_subset) - np.sum(shape1_subset) * np.sum(shape2_subset) / len(shape2_subset)
        denominator_1 = np.sum(np.square(shape1_subset)) - np.square(np.sum(shape1_subset)) / len(shape1_subset)
        denominator_2 = np.sum(np.square(shape2_subset)) - np.square(np.sum(shape2_subset)) / len(shape2_subset)
        denominator = 1
        if denominator_1 > 0 and denominator_2 > 0:
            denominator = np.sqrt(denominator_1 * denominator_2)
        else:
            numerator = 0
    else:
        raise ValueError("Maximum signal intensity not contained in region")
        
    #Return the cross-correlation result.
    return numerator / denominator
    
#Count signals above a value.
def count_above(threshold, annotation, signal, start, end, start_anno, end_anno, bin_sz):
    count = 0
    start_idx = start
    for sig in signal:
        is_between_anno = (start_anno <= start_idx) and (start_idx <= end_anno)
        if sig > threshold and (is_between_anno or annotation == ""):
            count += bin_sz
        start_idx += bin_sz
    return count