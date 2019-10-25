"""
Given the signal intensity and its association with ChromHMM, associate each region
with a regulatory element. The following arguments are required:
1. The input regions.
2. The list of signal associations with ChromHMM.
"""

import numpy as np
import scipy as sp
import sys
import os
import math
from scipy import stats
sys.path.append(os.path.abspath("../common_scripts"))
import pickle as pkl
import wig_and_signal_utils as wsu

def main():

    # Arguments
    input = sys.argv[1]
    signal = sys.argv[2]
    output = sys.argv[3]
    p_promoter = sys.argv[4]
    p_enhancer = sys.argv[5]
    p_weak = sys.argv[6]
    p_repressor = sys.argv[7]
   
    #Open all input files.
    regions = pkl.load(open(input, 'rb'))
    signals = pkl.load(open(signal, 'rb'))

    #Match the inputs to shapes and close the files.
    match_signals(regions, signals, output, p_promoter, p_enhancer, p_weak, p_repressor)
    
    """
Match each input with the given window size to the nearest shape for that window size.
Print out the region with its corresponding shape to a BED file.
"""
def match_signals(regions, signals, out, p_promoter, p_enhancer, p_weak, p_repressor):
    
    #Read in each line in the file and map it.
    out_f = open(out, "w")
    if(len(signals) > 1):
        for region in tqdm(regions):
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [label, confidence] = match_region(region, signals, p_promoter, p_enhancer, p_weak, p_repressor)
            
            #Build list of region matches.
            if label is not None:
                out_f.write("\t".join(region.chromosome, str(region.start), str(region.end), label, str(confidence)))
            
            
        #Print a message to the user.
        out_f.close()
        print("Files done")
    else:
        print("Only one shape. No annotation performed.")
        

"""
Find the closest match for the region in the list of signal
intensities. Do this by matching the median of the region
to its corresponding annotation.
"""
def match_region(region, signal_list, p_promoter, p_enhancer, p_weak, p_repressor):

    match = None
    
    #Find the intensity of the median signal within the region.
    # Use this to obtain the match.
    median = int(math.ceil(np.median(region.signals)))
    match = signal_list[median]
    
    # Get the annotation type comprising the maximum of signals
    # with this intensity.
    label = None
    argmax = wsu.get_annotation(match, p_promoter, p_enhancer, p_weak, p_repressor)
    if argmax == 0:
        label = "Promoter"
    elif argmax == 1:
        label = "Enhancer"
    elif argmax == 2:
        label = "Repressed"
    elif argmax == 3:
        label = "Weak"
        
    # Get the confidence of this annotation.
    total_anno = np.sum(match)
    percent_label = match[argmax] / total_anno

    return [label, percent_label]
 
    
if __name__ == "__main__":
    main()