"""
Given the shapes and their overlap with ChromHMM, associate each shape
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

BIN_SIZE = 50
def main():

    #Read in the bed file and shape data.
    bed = np.genfromtxt(sys.argv[1], delimiter='\t', dtype = str)
    shapes = pkl.load(open(sys.argv[2], "rb"))
    output = sys.argv[3]
    
    # Mapping from ChromHMM mnemonics to RE
    promoter = {"1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"}
    enhancer = {"6_EnhG", "7_Enh", "12_EnhBiv"}
    repressed = {"13_ReprPC", "14_ReprPCWk"}
    weak = {"9_Het", "15_Quies"}
        
    shape_col = 3 # Shape annotation
    bio_col = 8 # Biological (ChromHMM) annotation
    shape_start = 1 # Shape annotation start position
    shape_end = 2 # Shape annotation end position
    bio_start = 6 # ChromHMM annotation start position
    bio_end = 7 # ChromHMM annotation end position
    bio_len = 9 # ChromHMM annotation length
    

    #Get distribution of ChromHMM classes per shape.
    shape_names = get_names_of_shapes(shapes)
    total_percent_all = get_all_percentage_pairs(shape_col, bio_col, shape_start, shape_end, bio_start, bio_end,  bio_len, bed, promoter, enhancer, repressed, weak, shape_names)

    #Print all shapes with significant annotations, along with their annotations.
    save_associations(total_percent_all, shape, output)
    
#Compute percentage for each shape-annotation pair.
def get_all_percentage_pairs(anno, chrom_hmm_anno, start, end, chrom_hmm_start, chrom_hmm_end, chrom_hmm_len, bed, promoter, enhancer, repressed, weak, shapes):
    
    #Set up percentage matrix.
    sum_matrix = np.zeros((4, len(shapes)))
    cumulative_vec = np.zeros(len(shapes))
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    current_clust = "none"
    prev_start = -1
    prev_end = -1
    junk = 0
    
    for i in range(0, bed.shape[0]):
        
        #Get the previous data, if applicable.
        if i > 0:
            prev_start = int(bed[i - 1, start])
            prev_end = int(bed[i - 1, end])
            prev_clust = bed[i - 1, anno]
            
        #Get the next element data.
        next_line = bed[i,:]
        current_start = int(next_line[start])
        current_end = int(next_line[end])
        current_clust = next_line[anno]
        a = next_line[chrom_hmm_anno]
        idx = shapes.index(current_clust)
        
        #Add to the existing percentages.
        region = [float(i) for i in next_signal[3:len(next_signal)]]

        if a in promoter:
            sum_matrix[0, idx] += int(next_line[chrom_hmm_len])
        elif a in enhancer:
            sum_matrix[1, idx] += int(next_line[chrom_hmm_len])
        elif a in repressed:
            sum_matrix[2, idx] += int(next_line[chrom_hmm_len])
        elif a in weak:
            sum_matrix[3, idx] += int(next_line[chrom_hmm_len])
            
        #Add to the total sum for this shape.
        cumulative_vec[idx] = np.sum(sum_matrix[:,idx])
    
    #Get the set of percentages.
    cumulative_matrix = np.tile(cumulative_vec, (4, 1))
    return sum_matrix / cumulative_matrix
   
# For each shape, save its distribution.
def save_associations(percentages, shapes, out_name):
    
    #Get threshold to use in printing.
    out = open(out_name, 'a')
    perc_out = open(out_name + "_percents", 'a')
    percentage_out = open(out_name + "/" + chrom + "percentages", 'w')
    print("\n")
    scale = 5.0 / intensity
    perc_string = percentages.astype(str)
    
    #Print significant shapes with labels.
    labels = ["Promoter", "Enhancer", "Polycomb", "Weak"]
    for j in range(0, percentages.shape[1]):
        split_clust = shapes[j].split(",")
        scaled_clust = [float(i) for i in split_clust] * np.tile(scale, len(split_clust))
        scaled_clust_str = [str(i) for i in scaled_clust]
        perc_str = [str(i) for i in percentages[:,j]]
        joined = ','.join(scaled_clust_str)
        perc_joined = ','.join(perc_str)
        was_significant = False
        for i in range(0, percentages.shape[0]):
            other_percents = percentages[np.arange(percentages.shape[0])!=i,j]
            if labels[i] != "Polycomb" and percentages[i,j] > 0.5 and np.max(other_percents) <= percentages[i,j] / 2:
                out.write(labels[i] + "\t" + chrom + "_" + shape_names[j] + "\t" + joined)
                perc_out.write(labels[i] + "\t" + cell + "_" + chrom + "_" + shape_names[j] + "\t" + perc_joined)
                if joined.find('\n') != len(joined) - 1:
                    perc_out.write("\n")
                    out.write("\n")
                was_significant = True
        if not was_significant:
            out.write("Unknown" + "\t" + "_" + chrom + "_" + shape_names[j] + "\t" + joined)
            perc_out.write("Unknown" + "\t" + "_" + chrom + "_" + shape_names[j] + "\t" + perc_joined)
            if joined.find('\n') != len(joined) - 1:
                perc_out.write("\n")
                out.write("\n")
        #Print all percentages to a file to use later.
        percentage_out.write(','.join(perc_string[:,j]) + "\n")
            
    #Close the files.
    wig.close()
    out.close()
    percentage_out.close()
    
if __name__ == "__main__":
    main()