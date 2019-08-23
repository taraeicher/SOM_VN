"""
This script shifts the shapes learned from the SOM to merge shapes
separated only by phase. It outputs a file containing the filtered shapes.
"""
import sys
import os
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu
import region_defs
from tqdm import tqdm

"""
If one shape is a shift of another, merge them.
"""
def main():

    file_path = sys.argv[1]
    output_path = sys.argv[2]
    thresh = float(sys.argv[3])

    # Print message to user.
    print("Merging shifted shapes")

    # Open the file containing the SOM shapes.
    som_shapes = pkl.load(open(file_path, "rb"))
        
        #Compare each shape to all shapes after it.
        #If the shapes should be merged, shift the prior one as needed,
        #then delete the latter one.
        #A threshold of 0.75 is used like in CoSBI
        length = som_shapes.shape[0]
        prev_j = 0
        for i in tqdm(range(length - 1)):
            j = i + 1
            while j < length:
            
                #If i and j should be merged, merge them, get the new length,
                #and leave the comparison index as is. Else, increment.
                if should_merge(i, j, som_shapes, 0.75, thresh):
                    shift_shape(i, j, som_shapes)
                    length = len(som_shapes)
                else:
                    j += 1
                        
        #Save shifted shapes.
        pkl.dump(som_shapes, open(output_path, 'wb'))
    
#Determine whether two shapes should be merged.
def should_merge(idx1, idx2, shape_list, threshold, min_val):
    merge = False
    
    #Allow a shift of 1/4 of the shape dimension to each side.
    shift = 0.25 * len(shape_list[idx1].signals)
    
    #Get the maximum cross-correlation metric.
    max_cc = 0
    for delay in range(int(-shift), int(shift)):
        try:
            cross_correlation = wsu.get_crosscorr(shape_list[idx1].signals, shape_list[idx2].signals, delay, threshold, 2, True, True, min_val)
            if cross_correlation > max_cc:
                max_cc = cross_correlation

        #If the regions aren't valid for comparison at this delay, skip them.
        except ValueError:
            pass
            
    #If the maximum is above the threshold, then merge.

    if max_cc >= threshold:
        merge = True

    #Return decision about merging.
    return merge
    
#Use whichever shape has the most mass distributed at the center.
def shift_shape(shape, replacement, shape_list):
    
    #Determine which shape has its maximum closer to the center.
    max_index_shape = shape_list[shape].signals.index(max(shape_list[shape].signals))
    max_index_replacement = shape_list[replacement].signals.index(max(shape_list[replacement].signals)) 
    
    #If the latter shape has its max at the center, replace the prior shape with it.
    if abs((len(shape_list[replacement].signals) / 2) - max_index_replacement) < abs((len(shape_list[shape].signals) / 2) - max_index_shape):
        shape_list[shape] = shape_list[replacement]
        
    #Remove the latter shape.
    del shape_list[replacement]

if __name__ == "__main__":
    main()