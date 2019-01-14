import numpy as np
import scipy as sp
import sys
import os
import math
from scipy import stats
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu


"""
For each of the annotations, find its information gain for each shape.
For those shape-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the bed file and shape data.
    bed = np.genfromtxt(sys.argv[1], delimiter='\t', dtype = str)
    output = sys.argv[3]
    wig = sys.argv[4]
    chromosome = sys.argv[5]
    cell = sys.argv[6]
    anno_file = sys.argv[7]
    min_val = float(sys.argv[8])
        
    #shape annotations are in column 3. Biological annotations are in column 4.
    shape_col = 3
    bio_col = 8
    shape_start = 1
    shape_end = 2
    bio_start = 6
    bio_end = 7
    bio_len = 9
    
    #Get list of shapes.
    shapes = []
    shape_file = open(sys.argv[2], 'r')
    next_clust = shape_file.readline()
    while next_clust:
        shapes.append(next_clust)
        next_clust = shape_file.readline()
    unique_clusts = sorted(list(set(bed[:,shape_col])), key=lambda x: float(x))

    #Get percentage of each shape and each chromHMM annotation and each chromHMM annotation per shape.
    wig_file = open(wig, "r")
    threshold = wsu.get_intensity_percentile(0.75, wig_file, min_val)
    wig_file.close()
    total_percent_shapes = get_shape_percentages(shape_col, shape_start, shape_end, unique_clusts, bed)
    total_percent_anno = get_anno_percentages(bio_col, shape_start, shape_end, bio_len, bed)
    [total_percent_all, total_sums_all] = get_all_percentage_pairs(shape_col, bio_col, shape_start, shape_end, bio_start, bio_end,  bio_len, unique_clusts, bed, threshold, anno_file, shapes)

    #Print all shapes with significant annotations, along with their annotations.
    save_significant(total_percent_all, unique_clusts, shapes, wig, output, chromosome, cell, min_val)
   
#Get the percentage of the chromosome belonging to each shape.
def get_shape_percentages(anno, start, end, shapes, bed):

    #Set up dividends and divisor.
    sum_vec = np.zeros(len(shapes))
    cumulative = 0
    
    #Loop through bed file to compute percentage for each region.
    current_start = "-1"
    current_end = "-1"
    current_clust = "none"
    prev_start = "-1"
    prev_end = "-1"
    prev_clust = "none"
    
    for i in range(0, bed.shape[0]):
        
        #Get the previous data, if applicable.
        if i > 0:
            prev_start = bed[i - 1, start]
            prev_end = bed[i - 1, end]
            prev_clust = bed[i - 1, anno]
            
        #Get the next element data.
        next_line = bed[i,:]
        current_start = next_line[start]
        current_end = next_line[end]
        current_clust = next_line[anno]
        idx = shapes.index(current_clust)
            
        #Add to the total sum if the current start and end are not equal to the previous ones.
        if prev_start != current_start:
            sum_vec[idx] += int(current_end) - int(current_start)
            cumulative += int(current_end) - int(current_start)
    
    #Get the set of percentages.
    cumulative_vec = cumulative * np.ones(len(shapes))
    return sum_vec / cumulative_vec

#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def get_anno_percentages(chrom_hmm_anno, start, end, chrom_hmm_len, bed):

    #Set up percentage matrix.
    sum_vec = np.zeros(3)
    cumulative = 0
    
    #Loop through bed file to compute percentage for each region.
    current_start = "-1"
    current_end = "-1"
    current_clust = "none"
    prev_start = "-1"
    prev_end = "-1"
    
    for i in range(0, bed.shape[0]):
        
        #Get the previous data, if applicable.
        if i > 0:
            prev_start = bed[i - 1, start]
            prev_end = bed[i - 1, end]
            
        #Get the next element data.
        next_line = bed[i,:]
        current_start = next_line[start]
        current_end = next_line[end]
        a = next_line[chrom_hmm_anno]
        
        #Add to the existing percentages.
        if a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
            sum_vec[0] += int(next_line[chrom_hmm_len])
        elif next_line[chrom_hmm_len] != "0":
            sum_vec[1] += int(next_line[chrom_hmm_len])
        else:
            sum_vec[2] += int(current_end) - int(current_start)
            
        #Add to the total sum if the current start and end are not equal to the previous ones.
        if(prev_start != current_start and prev_start != "-1"):
            cumulative += int(prev_end) - int(prev_start)
    
    #Add the last cumulative sum.
    cumulative += int(prev_end) - int(prev_start)
    
    #Get the set of percentages.
    cumulative_vec = cumulative * np.ones(3)
    return sum_vec / cumulative_vec
    
#Compute percentage for each shape-annotation pair.
def get_all_percentage_pairs(anno, chrom_hmm_anno, start, end, chrom_hmm_start, chrom_hmm_end, chrom_hmm_len, shapes, bed, thresh, signals_path, shape_str):
    
    #Set up percentage matrix.
    sum_matrix = np.zeros((2, len(shapes)))
    cumulative_vec = np.zeros(len(shapes))
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    current_clust = "none"
    prev_start = -1
    prev_end = -1
    sigs = open(signals_path + "window3", "r")
    junk = 0
    
    for j in range(0,2):#Loop through BED file twice.
        next_signal = sigs.readline().split(",")
        if j == 1:
            next_signal = sigs.readline().split(",")
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

            clust_sig = [float(i) for i in shape_str[idx].split(",")]
            #Get the signal data.
            if (prev_start >= int(next_signal[1]) or current_start > int(next_signal[1])) and current_start > prev_start:
                next_signal = sigs.readline().split(",")
            if current_start > int(next_signal[1]):
                next_signal = sigs.readline().split(",")
            if int(next_signal[1]) == current_start:
                #Add to the existing percentages.
                region = [float(i) for i in next_signal[3:len(next_signal)]]
                count_clust = count_above(thresh, "", clust_sig, 0, len(clust_sig) * BIN_SIZE, 0, 0)
                count_a = count_above(thresh, a, region, current_start, current_end, int(next_line[chrom_hmm_start]), int(next_line[chrom_hmm_end]))
  
                if a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
                    sum_matrix[0, idx] += int(next_line[chrom_hmm_len]) if (count_clust == 0) else count_a
                elif a != "0":
                    sum_matrix[1, idx] += int(next_line[chrom_hmm_len]) if (count_clust == 0) else count_a
                
                #Add to the total sum if the current start and end are not equal to the previous ones.
                #if(prev_start != current_start and prev_start != "-1"):
            cumulative_vec[idx] = np.sum(sum_matrix[:,idx])
    
    #Get the set of percentages.
    cumulative_matrix = np.tile(cumulative_vec, (2, 1))
    return [sum_matrix / cumulative_matrix, np.sum(sum_matrix, 0)]
   
#Save the shapes that are significant.
def save_significant(percentages, shape_names, shapes, wig_name, out_name, chrom, cell, min_val):
    
    #Get threshold to use in printing.
    wig = open(wig_name, 'r')
    out = open(out_name, 'a')
    percentage_out = open(out_name + "_" + cell + "/" + chrom + "percentages", 'w')
    intensity = wsu.get_intensity_percentile(0.995, wig, min_val)
    print("\n")
    wig.close()
    scale = 5.0 / intensity
    perc_string = percentages.astype(str)
    
    #Print significant shapes with labels.
    labels = ["Enhancer", "Other"]
    for j in range(0, percentages.shape[1]):
        split_clust = shapes[j].split(",")
        scaled_clust = [float(i) for i in split_clust] * np.tile(scale, len(split_clust))
        scaled_clust_str = [str(i) for i in scaled_clust]
        joined = ','.join(scaled_clust_str)
        was_significant = False
        for i in range(0, percentages.shape[0]):
            other_percents = percentages[np.arange(percentages.shape[0])!=i,j]
            if percentages[i,j] > 0.5:
                out.write(labels[i] + "\t" + cell + "_" + chrom + "_" + shape_names[j] + "\t" + joined)
                if joined.find('\n') != len(joined) - 1:
                    out.write("\n")
                was_significant = True
        if not was_significant:
            out.write("Unknown" + "\t" + cell + "_" + chrom + "_" + shape_names[j] + "\t" + joined)
            if joined.find('\n') != len(joined) - 1:
                out.write("\n")
        #Print all percentages to a file to use later.
        percentage_out.write(','.join(perc_string[:,j]) + "\n")
            
    #Close the files.
    wig.close()
    out.close()
    percentage_out.close()

def count_above(threshold, annotation, signal, start, end, start_anno, end_anno):
    count = 0
    start_idx = start
    for sig in signal:
        is_between_anno = start_anno <= start_idx and start_idx <= end_anno
        if sig > threshold and (is_between_anno or annotation == ""):
            count += BIN_SIZE
        start_idx += BIN_SIZE
    return count
    
if __name__ == "__main__":
    main()