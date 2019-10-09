"""
Given the signal intensity and its overlap with ChromHMM, associate each intensity
with a regulatory element. Regulatory elements are defined as follows:
1. Promoter (Active, flanking active, bivalent, poised, and flanking
bivalent TSS)
2. Enhancer (Genic enhancers, enhancers, and bivalent enhancers)
3. Weak (Heterochromatin and quiescent)
4. Repressor (Repressed polycomb and weak repressed polycomb)
"""

import numpy as np
import scipy as sp
import sys
import os
import math
from scipy import stats
sys.path.append(os.path.abspath("../common_scripts"))
import pickle as pkl
from tqdm import tqdm

BIN_SIZE = 50
def main():

    signal_col = 4 # Signal intensity annotation
    bio_col = 8 # Biological (ChromHMM) annotation
    bin_start = 1 # Bin start position
    bin_end = 2 # Bin end position
    bio_len = 9 # Overlap length
    
    #Read in the bed file and shape data.
    bed = np.genfromtxt(sys.argv[1], delimiter='\t', dtype = str)
    signal_bins = get_all_bins(bed, signal_col)
    output = sys.argv[2]
    
    # Mapping from ChromHMM mnemonics to RE
    promoter = {"1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"}
    enhancer = {"6_EnhG", "7_Enh", "12_EnhBiv"}
    repressed = {"13_ReprPC", "14_ReprPCWk"}
    weak = {"9_Het", "15_Quies"}    

    #Get distribution of ChromHMM classes per signal bin.
    total_percent_all = get_all_percentage_pairs(signal_col, bio_col, bin_start, bin_end, bio_len, bed, promoter, enhancer, repressed, weak, signal_bins)
    pkl.dump(total_percent_all, open(output, "wb"))
    
"""
Extract all bins (at a resolution of 1 RPKM) from the BED file.   
Start with a large number of bins, then remove the unneeded bins. 
"""
def get_all_bins(bed, sig_c):

    max_threshold = 1000000
    bins = np.zeros(max_threshold)        
    bin_sz = 50
    #Use each entry in the file to calculate running metadata.
    print("Creating signal bins")
    for i in tqdm(range(bed.shape[0])):
        line = bed[i,:]
        val = float(line[sig_c])
        
        #Increment the appropriate location by the number of bins.
        bin = int(math.ceil(val))
        bins[bin] = bins[bin] + 1

    highest_bin = max_threshold - 1
    highest_bin_found = False
    bin = len(bins) - 1
    while bin > 0 and not highest_bin_found:
        if bins[bin] > 0:
            highest_bin_found = True
            highest_bin = bin
        else:
            bin = bin - 1

    return bins[0:highest_bin+1]
"""
Compute percentage for each shape-annotation pair.
"""
def get_all_percentage_pairs(anno, chrom_hmm_anno, start, end, chrom_hmm_len, bed, promoter, enhancer, repressed, weak, signal):
    
    #Set up percentage matrix.
    sum_matrix = np.zeros((4, len(signal)))
    cumulative_vec = np.zeros(len(signal))
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    current_shape = "none"
    prev_start = -1
    prev_end = -1
    print("Computing percentages")
    for i in tqdm(range(0, bed.shape[0])):
        
        #Get the previous data, if applicable.
        if i > 0:
            prev_start = int(bed[i - 1, start])
            prev_end = int(bed[i - 1, end])
            
        #Get the next element data.
        next_line = bed[i,:]
        current_start = int(next_line[start])
        current_end = int(next_line[end])
        current_signal = float(next_line[anno])
        chromhmm_annotation = next_line[chrom_hmm_anno]
        idx = math.ceil(current_signal)

        if chromhmm_annotation in promoter:
            sum_matrix[0, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation in enhancer:
            sum_matrix[1, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation in repressed:
            sum_matrix[2, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation in weak:
            sum_matrix[3, idx] += int(next_line[chrom_hmm_len])
    
    #Get the set of percentages.
    sum_totals = np.sum(sum_matrix, axis = 0)
    sum_totals_rep = np.clip(np.tile(sum_totals, (4,1)), a_min = 0.0001, a_max = None)
    return sum_matrix / sum_totals_rep
    
if __name__ == "__main__":
    main()