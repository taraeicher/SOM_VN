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

PROMOTER=0
ENHANCER=1
REPRESSOR=2
WEAK=3
OTHER=4

def main():

    # Arguments
    input = sys.argv[1]
    signal = sys.argv[2]
    output = sys.argv[3]
    is_peas = sys.argv[4] == "True"
   
    #Open all input files.
    regions = pkl.load(open(input, 'rb'))
    signals = pkl.load(open(signal, 'rb'))
    
    #Normalize signals by percentage of RE.
    #norm_signals = normalize_signals(signals)

    #Match the inputs to shapes and close the files.
    match_signals(regions, signals, output, is_peas)
   
"""
Normalize signal-RE associations by the sum of all signal associations for that RE,
thus obtaining the probability of the signal occurring given the RE. This is a maximum
likelihood distribution.   
"""
def normalize_signals(signals):
    normalized_signals = np.copy(signals)
    for i in range(signals.shape[0]):
        normalized_signals[i,:] = signals[i,:] / np.sum(signals, axis = 1)[i]
    np.nan_to_num(normalized_signals, copy=False)
    return normalized_signals
"""
Match each input with the given window size to the nearest shape for that window size.
Print out the region with its corresponding shape to a BED file.
"""
def match_signals(regions, signals, out, is_peas):
    
    #Read in each line in the file and map it.
    out_f = open(out, "w")
    if(len(signals) > 1):
        for region in tqdm(regions):
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [label, confidence] = match_region(region, signals, is_peas)
            
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
intensities. Do this by matching the max of the region
to its corresponding annotation.
"""
def match_region(region, signal_list, is_peas):

    match = None
    
    #Find the intensity of the median signal within the region.
    # Use this to obtain the match.
    argmax = int(math.ceil(np.amax(region.signals)))
    if argmax >= signal_list.shape[1]:
        argmax = signal_list.shape[1] - 1
    match = signal_list[:,argmax]
    labelval = 0
    
    # Get the annotation type comprising the maximum of signals
    # with this intensity.      
    label = None
    if is_peas:
        if match[ENHANCER] > match[OTHER]:
            label = "Enhancer"
            labelval = match[ENHANCER]
        else:
            label = "Other"
            labelval = match[OTHER]
    else:
        if np.argmax(match) == PROMOTER:
            label = "Promoter"
            labelval = match[PROMOTER]
        elif np.argmax(match) == ENHANCER:
            label = "Enhancer"
            labelval = match[ENHANCER]
        elif np.argmax(match) == REPRESSOR:
            label = "Repressor"
            labelval = match[REPRESSOR]
        else:
            label = "Weak"
            labelval = match[WEAK]

    # Get the confidence of this annotation.
    total_anno = np.sum(match)
    percent_label = labelval / total_anno

    return [label, percent_label]
 
    
if __name__ == "__main__":
    main()