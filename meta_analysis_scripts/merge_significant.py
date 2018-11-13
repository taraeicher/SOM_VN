from scipy import cluster
import numpy as np
import sys
import math
import common_ops as co
import os

"""
If one cluster is a shift of another, merge them.
"""
def main():

    file_path = sys.argv[1]
    output_path = sys.argv[2]
    merge_log = sys.argv[3]

    #Open the file containing the significant centroids.
    som_centroids = []
    som_anno = []
    som_names = []
    file = open(file_path, 'r')
    mlog = open(merge_log, 'w')
    next_line = file.readline()
    while next_line:
        split_tabs = next_line.split("\t")
        som_centroids.append([float(i) for i in split_tabs[2].split(",")])
        som_anno.append(split_tabs[0])
        som_names.append(split_tabs[1])
        next_line = file.readline()
            
    #Compare each cluster to all clusters after it.
    #If the clusters should be merged and have the same annotation, shift the prior one as needed,
    #then delete the latter one.
    #If the clusters should be merged and one is unknown, keep the one that is not unknown.
    #If the clusters should be merged and have different annotations, remove both.
    #A threshold of 0.75 is used like in CoSBI
    length = len(som_centroids)
    prev_j = 0
    for i in range(length - 1):
        j = i + 1
        while j < length:
            #If i and j meet the threshold, follow merging procedure.
            if should_merge(i, j, som_centroids, 0.75):
                mlog.write("Match between " + som_anno[i] + " and " + som_anno[j])
                mlog.write("(" + som_names[i] + ")" + " " + "(" + som_names[j] + ")")
                mlog.write("\n")
                #If they have the same annotation, merge them, get the new length,
                #and leave the comparison index as is.
                if som_anno[i] == som_anno[j]:
                    shift_cluster(i, j, som_centroids, som_names, som_anno)
                    length = len(som_centroids)
                #If one in unknown, keep the one that is not.
                #If they have different annotations, remove both to eliminate ambiguity.
                elif som_anno[i] == "Unknown":
                    del som_centroids[i]
                    del som_anno[i]
                    del som_names[i]
                    length = len(som_centroids)
                elif som_anno[j] == "Unknown":
                    del som_centroids[j]
                    del som_anno[j]
                    del som_names[j]
                    length = len(som_centroids)
                else:
                    del som_centroids[j]
                    del som_anno[j]
                    del som_names[j]
                    del som_centroids[i]
                    del som_anno[i]
                    del som_names[i]
                    length = len(som_centroids)
            else:
                j += 1
                        
        #Print shifted cluster centroids.
        outfile = open(output_path, 'w')
        for cluster in range(len(som_centroids)):
            outfile.write(som_names[cluster] + "\t" + som_anno[cluster] + "\t")
            for val in range(len(som_centroids[cluster])):
                outfile.write(str(som_centroids[cluster][val]))
                if val < len(som_centroids[cluster]) - 1:
                    outfile.write(",")
            outfile.write("\n")
        
    #Print message to user.
    print("Merging complete for all clusters in " + file_path)
    mlog.close()
    outfile.close()
    file.close()
    
#Determine whether two clusters should be merged.
def should_merge(clust1, clust2, cluster_list, threshold):
    merge = False

    #Allow a shift of 1/4 of the cluster dimension to each side.
    shift = 0.25 * len(cluster_list[clust1])
    #Get the maximum cross-correlation metric.
    max_cc = 0
    for delay in range(int(-shift), int(shift)):
        try:
            cross_correlation = co.get_crosscorr(cluster_list[clust1], cluster_list[clust2], delay, threshold, 2, True, True)
            if cross_correlation > max_cc:
                max_cc = cross_correlation

        #If the regions aren't valid for comparison at this delay, skip them.
        except ValueError:
            pass

    #If the maximum is above the threshold, then merge.
    if max_cc >= threshold:
        merge = True

    #Return decision about merging.
    return merge
    
#Use whichever cluster has the most mass distributed at the center.
def shift_cluster(cluster, replacement, cluster_list, names, annotations):
    
    #Determine which cluster has its maximum closer to the center.
    max_index_cluster = cluster_list[cluster].index(max(cluster_list[cluster]))
    max_index_replacement = cluster_list[replacement].index(max(cluster_list[replacement])) 
    
    #If the latter cluster has its max at the center, replace the prior cluster with it.
    if abs((len(cluster_list[replacement]) / 2) - max_index_replacement) < abs((len(cluster_list[cluster]) / 2) - max_index_cluster):
        cluster_list[cluster] = cluster_list[replacement]
        
    #Remove the latter cluster.
    del cluster_list[replacement]
    del annotations[replacement]
    del names[replacement]

if __name__ == "__main__":
    main()