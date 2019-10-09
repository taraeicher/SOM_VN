"""
Given the final set of learned shapes, annotate each signal region
with its best-matching shape. Output the annotation of each region
to a BED file. This script must be run on each chromosome.

The following arguments are required:
1. The input region file for the chromosome
2. The shape file for the chromosome
3. The output BED file for the chromosome.
"""

#Import sys for obtaining command line args.
import sys
import os
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs
import wig_and_signal_utils as wsu
import pandas as pd
import pickle as pkl
import numpy as np
from tqdm import tqdm

"""
Match each input to its closest shape. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    # Arguments
    input = sys.argv[1]
    shape = sys.argv[2]
    output = sys.argv[3]
   
    #Open all input files.
    regions = pkl.load(open(input, 'rb'))
    shapes = pkl.load(open(shape, 'rb'))

    #Match the inputs to shapes and close the files.
    match_shapes(regions, shapes, output)
    
"""
Match each input with the given window size to the nearest shape for that window size.
Print out the region with its corresponding shape to a BED file.
"""
def match_shapes(regions, shapes, out):
    
    #Read in each line in the file and map it.
    matched_regions = []
    if(len(shapes) > 1):
        for region in tqdm(regions):
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [match, crosscorr] = match_region(region, shapes)
            
            #Build list of region matches.
            if match is not None:
                matched_regions.append(region_defs.Matched_Region(region, shapes[match], crosscorr))
            
        save_to_bed(out, matched_regions)
            
        #Print a message to the user.
        print("Files done")
    else:
        print("Only one shape. No annotation performed.")
 
"""
Save a list of mapped regions to a BED file, where the ambiguity is the score of the match.
"""
def save_to_bed(bed, regions):
    bed_out = open(bed, "w")
    for region in regions:
        bed_out.write("chr" + str(region.chromosome) + "\t" + str(int(region.start)) + "\t" + str(int(region.end))
            + "\t" + str(region.shape.name) + "\t" + str(region.crosscorr) + "\n")
            
    bed_out.close()
 
"""
Find the closest match for the region in the list of shapes.
Return an ambiguity metric which measures how close the region
is to its nearest shape as opposed to other shapes.
"""
def match_region(region, shapes):

    #Create array to hold distances between region and shapes.
    max_crosscorr = 0
    match = None
    opt_delay = 0
    crosscorr_list = []
    
    #For each shape, determine the region's distance from it.
    for i in range(len(shapes)):
        shape = shapes[i].signals
        crosscorr, d = get_max_crosscorr(region.signals, shape)
        crosscorr_list.append(crosscorr)
        if crosscorr_list[i] > max_crosscorr:
            max_crosscorr = crosscorr_list[i]
            match = i
            opt_delay = d
            
    #The ambiguity metric is the ratio of the closest shape's distance
    #to the distance of the second closest shape.
    return [match, max_crosscorr]
    
"""
Find the minimum distance between the region and shifted shape.
"""
def get_max_crosscorr(region, shape):

    #Since the region is twice the size of the shape, search along the region and match at the best point.
    maximum = 0
    crosscorr = 0
    opt_delay = 0
    step = len(shape) / 10
    for i in range(0, len(shape)):
        try:
            delay = i * step
            crosscorr = wsu.get_crosscorr(region, shape, int(delay))
  
            #If the distance is the smallest so far, update.
            if crosscorr > maximum:
                maximum = crosscorr
                opt_delay = delay
        except ValueError:
            pass
    #Return the smallest distance
    return [crosscorr, delay]
    
if __name__ == "__main__":
    main()