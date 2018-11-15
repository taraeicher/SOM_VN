#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

#Import numpy for calculations
import numpy as np

#Import sys for obtaining command line args.
import sys

import common_ops as co

"""
For each datum in each of the window size files, match it to its
closest cluster in the cluster files. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    #Window size list and list of chromosomes are variable.
    #Allow user to choose number of iterations, initial learning
    #rate, and initial neighborhood size.
    #20, 50, 100, 250, 
    input = sys.argv[1]
    cluster = sys.argv[2]
    output = sys.argv[3]
    minimum = float(sys.argv[4])
    #clust_out_dir = sys.argv[4]

    #Create grids for labeling SOM nodes with cluster indices.
    som_centroids = []
    som_centroid_counts = []
    cluster_names = []
    cluster_count = 0
    
    #Open all input files.
    in_file = open(input, 'r')
    cluster_file = open(cluster, 'r')
    #cluster_out_files = [open(clust_out_dir + file.split(input_dir)[1], 'w') for file in sorted(glob.glob(input_dir + "*"))]

    #Notify the user that files were found
    print("Using the input file:")
    print(in_file.name)

    #Notify the user that files were found
    print("Using the cluster files:")
    print(cluster_file.name)

    #Match the inputs to clusters and close the files.
    match_clusters(in_file.name, cluster_file.name, output, minimum)
    in_file.close()
    cluster_file.close()
    
"""
Match each input with the given window size to the nearest cluster for that window size.
Print out the region with its corresponding cluster to a BED file.
"""
def match_clusters(in_file_name, cluster_file_name, out_dir, min_val):
    
    #Open input and cluster files.
    in_file = open(in_file_name, "r")
    cluster_file = open(cluster_file_name, "r")
    
    #Open output file.
    out_file = open(out_dir, "w")
    #cluster_file_out = open(cluster_out + "window" + str(win), "w")
    
    #Read in cluster data.
    clusters = []
    next_cluster = cluster_file.readline()
    counter = 0
    while next_cluster:
        clusters.append([float(i) for i in next_cluster.split(",")])
        next_cluster = cluster_file.readline()
    cluster_file.close()
    
    #Read in each line in the file and map it.
    if(len(clusters) > 1):
        next_line = in_file.readline()
        while(next_line):
        
            #Get the chromosome and position info to print to the output file.
            #Get the input signal to match with the clusters.
            inputStr = []
            labels = []
            split_line = next_line.split(",")
            labels = list(split_line[0:3])
            inputStr = list(split_line[3:len(split_line)])
            input = [float(i) for i in inputStr]
            
            #Match the data to the nearest cluster and obtain the match and the ambiguity metric.
            [match, ambig, out_str] = match_datum(input, clusters, min_val)
            
            #Print match to BED file. Format is:
            #chrom  start   end cluster_num 1 - ambiguity
            #Score is the opposite of the ambiguity metric.
            out_file.write("chr" + labels[0] + "\t" + labels[1] + "\t" + labels[2] + "\t" + str(match) + "\t" + str(1 - ambig) + "\t" + out_str + "\n")

            #Read the next line in the file.            
            next_line = in_file.readline()
            
        #Print a message to the user.
        print("Files done")
    else:
        print("Only one cluster. No annotation performed.")
    out_file.close()
    in_file.close()

"""
Find the closest match for the datum in the list of clusters.
Return an ambiguity metric which measures how close the datum
is to its nearest cluster as opposed to other clusters.
"""
def match_datum(datum, clusters, min_val):

    #Create array to hold distances between datum and clusters.
    max_crosscorr = 0
    match = 0
    opt_delay = 0
    crosscorr_list = []
    
    #For each cluster, determine the datum's distance from it.
    for i in range(len(clusters)):
        cluster = clusters[i]
        crosscorr, d = get_max_crosscorr(datum, cluster, min_val)
        crosscorr_list.append(crosscorr)
        if crosscorr_list[i] > max_crosscorr:
            max_crosscorr = crosscorr_list[i]
            match = i
            opt_delay = d

    #Calculate the ambiguity metric. Only do so if not all cross-correlations are zero.
    #If they are all zero, set ambiguity to 1 (completely ambiguous) and assign it to the
    #last cluster (artificially created, with all 0.1).
    ambig = 1.0
    if not max_crosscorr == 0:
        ambig = get_ambiguity(crosscorr_list, datum)
    else:
        match = len(clusters) - 1
    
    #Write the cluster values to a separate file.
    #This file will be used for cluster validity analysis.
    
    comma = ","
    out_str = comma.join(str(e) for e in datum[int(opt_delay):len(clusters[0]) + int(opt_delay)])
    #cluster_out_file.write(out_str + "\n")
            
    #The ambiguity metric is the ratio of the closest cluster's distance
    #to the distance of the second closest cluster.
    return [match, ambig, out_str]

"""
Compute the ambiguity metric. It is equivalent to the smallest distance
from the datum to a cluster centroid, divided by the second smallest
distance from the datum to a cluster centroid.
""" 
def get_ambiguity(list, datum):
    list_array = np.asarray(list)
    max_idx = np.argmax(list_array)
    list_max_removed = np.delete(list_array, max_idx)
    max_idx_2 = np.argmax(list_max_removed)
    ambiguity = list_max_removed[max_idx_2] / list_array[max_idx]
    return ambiguity

"""
Find the minimum distance between the datum and shifted cluster.
"""
def get_max_crosscorr(datum, cluster, min_val):

    #Since the datum is twice the size of the cluster, search along the datum and match at the best point.
    maximum = 0
    crosscorr = 0
    opt_delay = 0
    step = len(cluster) / 10
    for i in range(0, 10):
        try:
            delay = i * step
            crosscorr = co.get_crosscorr(datum, cluster, int(delay), 0.75, 0, False, False, min_val)
                
            #If the distance is the smallest so far, update.
            if crosscorr > maximum:
                maximum = crosscorr
                opt_delay = delay
        except ValueError:
            pass
    #Return the smallest distance
    return [crosscorr, delay]
    
"""
Find the minimum distance between the datum and shifted cluster.
"""
def get_min_shifted_distance(datum, cluster):

    #Since the datum is twice the size of the cluster, search along the datum and match at the best point.
    min_dist = sum(datum)
    for delay in range(0, len(cluster)):
        datum_array = np.asarray(datum[delay:len(cluster) + delay])
        clust_array = np.asarray(cluster)
        dist = np.linalg.norm(clust_array - datum_array)
            
        #If the distance is the smallest so far, update.
        if dist < min_dist and np.max(np.asarray(datum)) == np.max(datum_array):
            min_dist = dist
    #Return the smallest distance
    return min_dist

if __name__ == "__main__":
    main()