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
from tqdm import tqdm
sys.path.append(os.path.abspath("../common_scripts"))
import pickle as pkl
import wig_and_signal_utils as wsu

def main():

    # Arguments
    input = sys.argv[1]
    signal = sys.argv[2]
    output = sys.argv[3]
    p_promoter = float(sys.argv[4])
    p_enhancer = float(sys.argv[5])
    p_repressor = float(sys.argv[6])
    p_weak = float(sys.argv[7])
    is_peas = sys.argv[8] == "True"
   
    #Open all input files.
    regions = pkl.load(open(input, 'rb'))
    signals = pkl.load(open(signal, 'rb'))

    #Match the inputs to shapes and close the files.
    match_signals(regions, signals, output, p_promoter, p_enhancer, p_weak, p_repressor, is_peas)
    
    """
Match each input with the given window size to the nearest shape for that window size.
Print out the region with its corresponding shape to a BED file.
"""
def match_signals(regions, signals, out, p_promoter, p_enhancer, p_weak, p_repressor, is_peas):
    
    #Read in each line in the file and map it.
    out_f = open(out, "w")
    if(len(signals) > 1):
        for region in tqdm(regions):
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [label, confidence] = match_region(region, signals, p_promoter, p_enhancer, p_weak, p_repressor, is_peas)
            
            #Build list of region matches.
            if label is not None:
                out_f.write("\t".join(["chr" + region.chromosome, str(int(region.start)), str(int(region.end)), label, str(confidence)]) + "\n")
            
            
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
def match_region(region, signal_list, p_promoter, p_enhancer, p_weak, p_repressor, is_peas):

    match = None
    
    #Find the intensity of the median signal within the region.
    # Use this to obtain the match.
    argmax = int(math.ceil(np.amax(region.signals)))#int(math.ceil(np.median(region.signals)))
    match = signal_list[:,argmax]
    labelval = 0
    
    # Get the annotation type comprising the maximum of signals
    # with this intensity.      
    label = None
    if match[0] >= p_promoter and match[1] < p_enhancer:
        label = "Promoter"
        labelval = match[0]
    elif match[1] >= p_enhancer and match[0] < p_promoter:
        label = "Enhancer"
        labelval = match[1]
    elif match[1] >= p_enhancer and match[0] >= p_promoter:
        if match[1] >= match[0]:
            label = "Enhancer"
            labelval = match[1]
        else:
            label = "Promoter"
            labelval = match[0]
    elif is_peas:
        label = "Other"
        labelval = match[4]
    elif match[2] >= match[3]:
        label = "Repressor"
        labelval = match[2]
    else:
        label = "Weak"
        labelval = match[4]
        
    # Get the confidence of this annotation.
    total_anno = np.sum(match)
    percent_label = labelval / total_anno

    return [label, percent_label]
 
    
if __name__ == "__main__":
    main()