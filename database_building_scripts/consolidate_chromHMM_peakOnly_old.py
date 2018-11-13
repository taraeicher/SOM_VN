import numpy as np
import scipy as sp
import sys
import os
import math
import common_ops as ops
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc


"""
For each of the annotations, find its information gain for each cluster.
For those cluster-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the bed file and cluster data.
    bed = np.genfromtxt(sys.argv[1], delimiter='\t', dtype = str)
    output = sys.argv[3]
    wig = sys.argv[4]
    img_dir = sys.argv[5]
    chromosome = sys.argv[6]
    window_size = sys.argv[7]
    out_pvals = sys.argv[8]
    anno_file = sys.argv[9]
    cell = sys.argv[10]
        
    #Cluster annotations are in column 3. Biological annotations are in column 4.
    cluster_col = 3
    bio_col = 8
    cluster_start = 1
    cluster_end = 2
    bio_start = 6
    bio_end = 7
    bio_len = 9
    
    #Get list of clusters.
    clusters = []
    cluster_file = open(sys.argv[2], 'r')
    next_clust = cluster_file.readline()
    while next_clust:
        clusters.append(next_clust)
        next_clust = cluster_file.readline()
    unique_clusts = sorted(list(set(bed[:,cluster_col])), key=lambda x: float(x))

    #Get percentage of each cluster and each chromHMM annotation and each chromHMM annotation per cluster.
    wig_file = open(wig, "r")
    threshold = ops.get_intensity_percentile(0.75, wig_file)
    wig_file.close()
    total_percent_clusters = get_cluster_percentages(cluster_col, cluster_start, cluster_end, unique_clusts, bed)
    total_percent_anno = get_anno_percentages(bio_col, cluster_start, cluster_end, bio_len, bed)
    [total_percent_all, total_sums_all] = get_all_percentage_pairs(cluster_col, bio_col, cluster_start, cluster_end, bio_start, bio_end,  bio_len, unique_clusts, bed, threshold, anno_file, clusters)

    #Print all clusters with significant annotations, along with their annotations.
    #print_cluster_stats(total_percent_all, total_sums_all, total_percent_anno, unique_clusts, clusters, wig, output, out_pvals)
    save_significant(total_percent_all, unique_clusts, clusters, wig, output, chromosome, cell)
    #save_line_charts(total_percent_all, total_percent_anno, unique_clusts, img_dir, chromosome, window_size)
   
#Get the percentage of the chromosome belonging to each cluster.
def get_cluster_percentages(anno, start, end, clusters, bed):

    #Set up dividends and divisor.
    sum_vec = np.zeros(len(clusters))
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
        idx = clusters.index(current_clust)
            
        #Add to the total sum if the current start and end are not equal to the previous ones.
        if prev_start != current_start:
            sum_vec[idx] += int(current_end) - int(current_start)
            cumulative += int(current_end) - int(current_start)
    
    #Get the set of percentages.
    cumulative_vec = cumulative * np.ones(len(clusters))
    return sum_vec / cumulative_vec

#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def get_anno_percentages(chrom_hmm_anno, start, end, chrom_hmm_len, bed):

    #Set up percentage matrix.
    sum_vec = np.zeros(5)
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
        if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
            sum_vec[0] += int(next_line[chrom_hmm_len])
        elif a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
            sum_vec[1] += int(next_line[chrom_hmm_len])
        elif a == "13_ReprPC" or a == "ReprPCWk":
            sum_vec[2] += int(next_line[chrom_hmm_len])
        elif a == "9_Het" or a == "15_Quies":
            sum_vec[3] += int(next_line[chrom_hmm_len])
        elif next_line[chrom_hmm_len] != "0":
            sum_vec[4] += int(next_line[chrom_hmm_len])
        else:
            sum_vec[4] += int(current_end) - int(current_start)
            
        #Add to the total sum if the current start and end are not equal to the previous ones.
        if(prev_start != current_start and prev_start != "-1"):
            cumulative += int(prev_end) - int(prev_start)
    
    #Add the last cumulative sum.
    cumulative += int(prev_end) - int(prev_start)
    
    #Get the set of percentages.
    cumulative_vec = cumulative * np.ones(5)
    return sum_vec / cumulative_vec
    
#Compute percentage for each cluster-annotation pair.
def get_all_percentage_pairs(anno, chrom_hmm_anno, start, end, chrom_hmm_start, chrom_hmm_end, chrom_hmm_len, clusters, bed, thresh, signals_path, cluster_str):
    
    #Set up percentage matrix.
    sum_matrix = np.zeros((4, len(clusters)))
    cumulative_vec = np.zeros(len(clusters))
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    current_clust = "none"
    prev_start = -1
    prev_end = -1
    sigs = open(signals_path, "r")
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
            idx = clusters.index(current_clust)

            #print(current_clust + " " + str(len(cluster_str)))
            clust_sig = [float(i) for i in cluster_str[idx].split(",")]
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
  
                if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
                    sum_matrix[0, idx] += int(next_line[chrom_hmm_len]) if (count_clust == 0) else count_a
                elif a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
                    sum_matrix[1, idx] += int(next_line[chrom_hmm_len]) if (count_clust == 0) else count_a
                elif a == "13_ReprPC" or a == "ReprPCWk":
                    sum_matrix[2, idx] += int(next_line[chrom_hmm_len]) if (count_clust == 0) else count_a
                elif a == "9_Het" or a == "15_Quies":
                    sum_matrix[3, idx] += int(next_line[chrom_hmm_len]) if (count_clust == 0) else count_a
                #Add to the total sum if the current start and end are not equal to the previous ones.
                #if(prev_start != current_start and prev_start != "-1"):
            cumulative_vec[idx] = np.sum(sum_matrix[:,idx])
        
    #Add the last cumulative sum.
    #cumulative_vec[clusters.index(prev_clust)] += int(int(prev_end) - int(prev_start)) if (count_null == 0) else count_null
    
    #Get the set of percentages.
    cumulative_matrix = np.tile(cumulative_vec, (4, 1))
    return [sum_matrix / cumulative_matrix, np.sum(sum_matrix, 0)]
   
#Save the clusters that are significant.
def save_significant(percentages, cluster_names, clusters, wig_name, out_name, chrom, cell):
    
    #Get threshold to use in printing.
    wig = open(wig_name, 'r')
    out = open(out_name, 'a')
    percentage_out = open(out_name + "_" + cell + "/" + chrom + "percentages", 'w')
    intensity = ops.get_intensity_percentile(0.995, wig)
    print("\n")
    wig.close()
    scale = 5.0 / intensity
    perc_string = percentages.astype(str)
    
    #Print significant clusters with labels.
    labels = ["Promoter", "Enhancer", "Polycomb", "Weak"]
    for j in range(0, percentages.shape[1]):
        split_clust = clusters[j].split(",")
        scaled_clust = [float(i) for i in split_clust] * np.tile(scale, len(split_clust))
        scaled_clust_str = [str(i) for i in scaled_clust]
        joined = ','.join(scaled_clust_str)
        was_significant = False
        for i in range(0, percentages.shape[0]):
            other_percents = percentages[np.arange(percentages.shape[0])!=i,j]
            if percentages[i,j] > 0.5 and np.max(other_percents) <= percentages[i,j] / 2:
                out.write(labels[i] + "\t" + cell + "_" + chrom + "_" + cluster_names[j] + "\t" + joined)
                if joined.find('\n') != len(joined) - 1:
                    out.write("\n")
                was_significant = True
        if not was_significant:
            out.write("Unknown" + "\t" + cell + "_" + chrom + "_" + cluster_names[j] + "\t" + joined)
            if joined.find('\n') != len(joined) - 1:
                out.write("\n")
        #Print all percentages to a file to use later.
        percentage_out.write(','.join(perc_string[:,j]) + "\n")
            
    #Close the files.
    wig.close()
    out.close()
    percentage_out.close()

#Save a bar chart with the ChromHMM distribution for each cluster, compared to the baseline.  
def save_line_charts(percentages, anno_percentages, cluster_names, path, chr, w):

    # plot
    names = cluster_names
    greenline = np.asarray(list(percentages[0,:]))
    orangeline = np.asarray(list(percentages[1,:]))
    blueline = np.asarray(list(percentages[2,:]))
    redline = np.asarray(list(percentages[3,:]))
    r = range(0, len(names))
    plt.figure(figsize=(12,6)) 
    # Create green Bars
    plt.plot(r, greenline, color='#008000', label = "Promoter")
    # Create orange Bars
    plt.plot(r,  orangeline, color='#ffa500', label = "Enhancer")
    # Create blue Bars
    plt.plot(r, blueline, color='#0000ff', label = "Polycomb")
    # Create red Bars
    plt.plot(r, redline, color='#ff0000', label = "Weak")
    
    # Add title, axis, and legend. Save plot.
    plt.xticks(r, names)
    plt.legend()
    plt.xlabel("Shape")
    plt.ylabel("Percentage")
    plt.title("Chromosome " + chr)
    plt.savefig(path + chr + ".png")

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