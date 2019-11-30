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
import region_defs

def main():

    #Read in the bed file and shape data.
    bed = np.genfromtxt(sys.argv[1], delimiter='\t', dtype = str)
    shapes = pkl.load(open(sys.argv[2], "rb"))
    output = sys.argv[3]
    is_peas = sys.argv[4] == "True"
    
    # Mapping from ChromHMM mnemonics to RE
    promoter = {"1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"}
    enhancer = {"6_EnhG", "7_Enh", "12_EnhBiv"}
    repressed = {"13_ReprPC", "14_ReprPCWk"}
    weak = {"9_Het", "15_Quies"}
    
    if is_peas:
        promoter = {}
        enhancer = {"AE", "OE"}
        repressed = {}
        weak = {}
        
    shape_col = 3 # Shape annotation
    bio_col = 8 # Biological (ChromHMM) annotation
    shape_start = 1 # Shape annotation start position
    shape_end = 2 # Shape annotation end position
    bio_start = 6 # ChromHMM annotation start position
    bio_end = 7 # ChromHMM annotation end position
    bio_len = 9 # ChromHMM annotation length
    

    #Get distribution of ChromHMM classes per shape.
    shape_names = get_names_of_shapes(shapes)
    total_percent_all = get_all_percentage_pairs(shape_col, bio_col, shape_start, shape_end, bio_start, bio_end,  bio_len, bed, promoter, enhancer, repressed, weak, shape_names, is_peas)

    #Build the list of shapes.
    assoc_shapes = []
    for i in range(len(shapes)):
        shape = shapes[i]
        assoc_shapes.append(region_defs.Shape_Association(shape, total_percent_all[0,i], total_percent_all[1,i], total_percent_all[2,i], total_percent_all[3,i], total_percent_all[4,i]))
    
    #Print all shapes with significant annotations, along with their annotations.
    pkl.dump(assoc_shapes, open(output, "wb"))
    
"""
Compute percentage for each shape-annotation pair.
"""
def get_all_percentage_pairs(anno, chrom_hmm_anno, start, end, chrom_hmm_start, chrom_hmm_end, chrom_hmm_len, bed, promoter, enhancer, repressed, weak, shapes, is_peas):
    
    #Set up percentage matrix.
    sum_matrix = np.zeros((5, len(shapes)))
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    current_shape = "none"
    prev_start = -1
    prev_end = -1
    
    for i in range(0, bed.shape[0]):
        
        #Get the previous data, if applicable.
        if i > 0:
            prev_start = int(bed[i - 1, start])
            prev_end = int(bed[i - 1, end])
            
        #Get the next element data.
        next_line = bed[i,:]
        current_start = int(next_line[start])
        current_end = int(next_line[end])
        current_shape = next_line[anno]
        chromhmm_annotation = next_line[chrom_hmm_anno]
        idx = shapes.index(current_shape)

        if chromhmm_annotation in promoter:
            sum_matrix[0, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation in enhancer:
            sum_matrix[1, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation in repressed:
            sum_matrix[2, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation in weak:
            sum_matrix[3, idx] += int(next_line[chrom_hmm_len])
        elif chromhmm_annotation != "." and is_peas:
            sum_matrix[4, idx] += int(next_line[chrom_hmm_len])
    
    #Get the set of percentages.
    sum_totals = np.sum(sum_matrix, axis = 0)
    sum_totals_rep = np.clip(np.tile(sum_totals, (5,1)), a_min = 0.0001, a_max = None)
    return sum_matrix / sum_totals_rep
   
"""
Compile the names of all shapes.
"""
def get_names_of_shapes(shapes):
    shape_names = []
    for shape in shapes:
        shape_names.append(shape.name)
    return shape_names
    
if __name__ == "__main__":
    main()