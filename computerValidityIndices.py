#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

#Import numpy for calculations
import numpy as np
import scipy as sp

#Import glob for obtaining the files.
import glob, os

#Import sys for obtaining command line args.
import sys

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
	num_window_sizes = 7
	input_dir = sys.argv[1]
	
	#Open all input files.
	in_file_names = [[]]
	for window in range(0, num_window_sizes):
		in_file_names.add([input_dir + str(window + file for file in sorted(os.listdir(input_dir + window))])
	
		#Notify the user that files were found
		print("For window size " + str(window) + ", using the input files:")
		for file in in_file_names:
			print(file)

	#Find and print indices.
	Parallel(n_jobs=num_window_sizes)(delayed(find_indices)(in_file_names[i], i) for i in range(0, num_window_sizes))
	for file in in_files:
		file.close()
	
"""
Calculate and print the indices for each clustering metric.
"""
def find_indices(in_file_names, window):
	
	#Open input and cluster files.
	in_files = [open(name, "r") for name in in_file_names]
	
	#Add all clusters to a list, then create a list of regions for each cluster.
	#Structure is a list of lists.
	clusters = [[]]
	fill_in_data(clusters)
	cluster_data = np.array(np.array(np.array(clusters)))
	
	#Get components needed for indices.
	cluster_means = get_cluster_means(cluster_data)
	inter_cluster_distances = get_inter_cluster_mean_distances(cluster_means)
	total_mean = get_total_mean(cluster_data)
	cluster_variances = get_dimensionwise_cluster_variances(cluster_data, total_mean)
	cluster_variances_dimensionless = get_dimensionless_cluster_variances(cluster_data)
	dataset_variances = get_dimensionwise_dataset_variance(cluster_data, total_mean)
	cluster_counts = get_cluster_counts(cluster_data)
	variance = get_dataset_variance(cluster_data)
	ssw = get_ssw(cluster_data)
	
	#Calculate and print each index.
	print_dunn(cluster_data, window)
	print_davies_bouldin(cluster_means, cluster_counts, cluster_variances, inter_cluster_distances, window)
	print_rmsstd(len(cluster_data[0][0]), cluster_counts, ssw, window)
	print_rs(cluster_variances, ssw, dataset_variance)
	print_sd(cluster_variances, datset_variance, len(cluster_data[0]), inter_cluster_distances, window)
	print_s_dbw(cluster_data, window)
	
###################### CLUSTER VALIDITY INDICES ############################

"""
Print the Dunn Index result
"""
def print_dunn(cluster_data, window):

	#Calculate Dunn index.
	inter_cluster_distances = get_inter_cluster_min_distances(cluster_data)
	max_cluster_diam = np.max(get_max_cluster_diameter(cluster_data))
	dunn_index = np.min(get_min_pairwise_dist(inter_cluster_distances / max_cluster_diam), axis = 1)
	
	#Print the value.
	print("Dunn Index for window " + str(window) + " is " + str(dunn_index))

"""
Print the Davies-Bouldin Index.
"""
def print_davies_bouldin(cluster_means, cluster_counts, cluster_variances, mean_dists, window):

	#Calculate Davies-Bouldin index
	
	similarity_measures = get_similarity_measures(inter_cluster_dists, cluster_variances, cluster_counts)
	davies_bouldin_index = np.sum(np.max(similarity_measures, axis = 0)) / len(cluster_means)
	
	#Print the value.
	print("Davies-Bouldin Index for window " + str(window) + " is " + str(davies_bouldin_index))
	
"""
Print the RMSSTD index.
"""
def print_rmsstd(num_dimensions, cluster_counts, ssw):

	#Calculate RMSSTD index
	denominator = np.sum(cluster_counts * num_dimensions)
	rmsstd = np.sqrt(np.div(ssw, denominator))
	
	#Print the value.
	print("RMSSTD Index for window " + str(window) + " is " + str(rmsstd))
	
"""
Print the RS index.
"""
def print_rs(ssw, cluster_variances_dimensionless, window):

	#Calculate RS index.
	sst = np.sum(np.sum(cluster_variances_dimensionless, axis = 0), axis = 1)
	rs = (sst - ssw) / sst
	
	#Print the value.
	print("RS Index for window " + str(window) + " is " + str(rs))
	
"""
Print the SD index.
"""
def print_sd(cluster_variances, dataset_variances, num_clusters, window):

	#Calculate the SD index.
	scattering = get_scattering(cluster_variances, dataset_variances, num_clusters)
	total_separation = get_dis(inter_cluster_distances)
	alpha = 1
	sd = alpha * scattering + total_separation
	
	#Print the value.
	print("SD Index for window " + str(window) + " is " + str(sd))
	
############### TOOLS FOR CALCULATING INDICES #########################
"""
Return the sum of variances across dimensions and clusters.
"""
def get_ssw(cluster_variances)
	return np.sum(np.sum(cluster_variances, axis = 2))

"""
Return the mean of each cluster in each dimension.
"""
def get_cluster_means(cluster_data)
	return np.mean(cluster_data, axis = 1)
	
"""
Return the mean of each dimension across clusters.
"""
def get_dimensionwise_means(cluster_data)
	return np.mean(cluster_data, axis = (0, 1))
	
"""
Return the variance of each dimension across clusters.
"""
def get_dimensionwise_cluster_variances(cluster_data)
	return np.var(cluster_data, axis = (0, 1))
	
"""
Return the variance across clusters over all dimensions.
"""
def get_dimensionless_cluster_variances(cluster_data)
	return np.var(cluster_data, axis = 0)
	
"""
Return the variance of each dimension across the entire dataset.
"""
def get_dimensionwise_cluster_variances(cluster_data)
	return np.var(cluster_data, axis = 1)

"""
Return the number of data mapping to each cluster.
"""
def get_cluster_counts(cluster_data)
	return dim(cluster_data)[0]

"""
Find the minimum pairwise distance between elements for each pair of clusters.
"""
def get_inter_cluster_min_distances(cluster_data):
	
	#Get the minimum pairwise distance between elements in each pair of clusters.
	distances = np.min([get_min_pairwise_dist(cluster_data,i) for i in enumerate(cluster_data)])
	
"""
Get minimum distance between all pairs of values.
"""
def get_min_pairwise_dist(data, index)
	return np.min(sp.cdist(data[index], data[np.arange(len(data))!= index], euclidean))

"""
Get maximum distance over all pairs of elements in a cluster.
"""
def get_max_distance_in_cluster(data)
	return np.max(sp.cdist(data, data, euclidean))

"""
Return pairwise distances between all cluster means.
"""
def get_inter_cluster_mean_distances(cluster_means)
	return sp.cdist(cluster_means, cluster_means)

"""
Find the maximum pairwise distance between elements within each cluster.
"""
def get_max_cluster_diameter(cluster_data):
	max_diameter = np.apply_along_axis(get_max_distance_in_cluster, 0, cluster_data, cluster_data)

"""
Calculate similarity metrics in the following manner:
1. Calculate average distance from each cluster element in cluster i to cluster center of i.
2. For each pair i,j, sum together dispersions for clusters i and j and divide by distance between i and j.
"""
def get_similarity_measures(inter_clust_dists, cluster_variances, cluster_counts):
	
	#Calculate dispersion for each cluster. This is the standard deviation of the cluster divided by its size.
	dispersions = np.sqrt(cluster_variances) / cluster_counts
	
	#For each pair of clusters, return similarity measure.
	pairwise_dispersion_sums = np.apply_along_axis(get_sum_to_all, 0, dispersions, dispersions)
	
	#Return dispersions divided by pairwise distances.
	return pairwise_dispersion_sums / inter_clust_dists
	
"""
For each distance value, get the sum between it and all other distance values.
"""
def get_sum_to_all(value, vals)
	value_tile = np.tile(value, len(vals))
	return value_tile + vals
	
"""
Return the variance of each dimension across clusters.
"""
def get_scattering(cluster_variances, dataset_variances, num_clusters)
	return np.sum(np.linalg.norm(cluster_variances) / np.linalg.norm(dataset_variances)) / num_clusters

"""
Return the variance of each dimension across clusters.
"""
def get_dis(inter_clust_dists)
	sum_1_dim = 1 / np.sum(np.abs(inter_clust_dists), axis = 0)
	sum = np.sum(sum_1_dim)
	max_dist = np.max(np.abs(inter_clust_dists))
	min_dist = np.min(np.abs(inter_clust_dists)) #Find some way to remove diagonal.
	return (max_dist / min_dist) * sum
	
if __name__ == "__main__":
	main()