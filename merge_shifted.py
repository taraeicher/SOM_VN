from scipy import cluster
import numpy as np
import sys
import math
import common_ops as co

"""
If one cluster is a shift of another, merge them.
"""
def main():

	num_windows = int(sys.argv[1])
	file_path = sys.argv[2]
	output_path = sys.argv[3]

	#For each file, obtain the SOM centroids and cluster them.
	for win in range(0, num_windows):
	
		#Print message to user.
		print("Merging shifted clusters for window " + str(win))
	
		#Open the file containing the SOM centroids.
		#Add all centroids to the list of data to cluster.
		som_centroids = []
		file = open(file_path + str(win), 'r')
		next_line = file.readline()
		while next_line:
			som_centroids.append([float(i) for i in next_line.split(",")])
			next_line = file.readline()
		
		#Compare each cluster to all clusters after it.
		#If the clusters should be merged, shift the prior one as needed,
		#then delete the latter one.
		#A threshold of 0.75 is used like in CoSBI
		length = len(som_centroids)
		prev_j = 0
		for i in range(length - 1):
			j = i + 1
			while j < length:
			
				#If i and j should be merged, merge them, get the new length,
				#and leave the comparison index as is. Else, increment.
				if should_merge(i, j, som_centroids, 0.75):
					shift_cluster(i, j, som_centroids)
					length = len(som_centroids)
				else:
					j += 1
						
		#Print shifted cluster centroids.
		file = open(output_path + str(win), 'w')
		for cluster in range(len(som_centroids)):
			for val in range(len(som_centroids[cluster])):
				file.write(str(som_centroids[cluster][val]))
				if val < len(som_centroids[cluster]) - 1:
					file.write(",")
			file.write("\n")
		
	#Print message to user.
	print("Merging complete for all windows.")
	
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
def shift_cluster(cluster, replacement, cluster_list):
	
	#Determine which cluster has its maximum closer to the center.
	max_index_cluster = cluster_list[cluster].index(max(cluster_list[cluster]))
	max_index_replacement = cluster_list[replacement].index(max(cluster_list[replacement]))	
	
	#If the latter cluster has its max at the center, replace the prior cluster with it.
	if abs((len(cluster_list[replacement]) / 2) - max_index_replacement) < abs((len(cluster_list[cluster]) / 2) - max_index_cluster):
		cluster_list[cluster] = cluster_list[replacement]
		
	#Remove the latter cluster.
	del cluster_list[replacement]

if __name__ == "__main__":
	main()