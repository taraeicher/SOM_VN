#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

#Import numpy for calculations
import numpy as np

#Import glob for obtaining the files.
import glob, os

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
	num_window_sizes = 7
	input_dir = sys.argv[1]
	cluster_dir = sys.argv[2]
	output_dir = sys.argv[3]

	#Create grids for labeling SOM nodes with cluster indices.
	som_centroids = []
	som_centroid_counts = []
	for i in range(num_window_sizes):
		som_centroids.append(list())
		som_centroid_counts.append(list())
	cluster_names = [[]]
	cluster_count = 0
	
	#Open all input files.
	in_files = [open(input_dir + file, 'r') for file in sorted(os.listdir(input_dir))]
	cluster_files = [open(cluster_dir + file, 'r') for file in sorted(os.listdir(cluster_dir))]
	
	if len(in_files) == num_window_sizes:
		#Notify the user that files were found
		print("Using the input files:")
		for file in in_files:
			print(file.name)
			
		if len(cluster_files) == num_window_sizes:
			#Notify the user that files were found
			print("Using the cluster files:")
			for cfile in cluster_files:
				print(cfile.name)

			#Match the inputs to clusters and close the files.
			#Parallel(n_jobs=num_window_sizes)(delayed(match_clusters)(in_files[i].name, i, cluster_files[i].name, output_dir) for i in range(0, num_window_sizes))
			Parallel(n_jobs=2)(delayed(match_clusters)(in_files[i].name, i, cluster_files[i].name, output_dir) for i in range(0, num_window_sizes))
			for file in in_files:
				file.close()
			for file in cluster_files:
				cfile.close()
			
		#Print an error if number of files incorrect.
		else:
			print("Error: number of files in cluster directory not equal to number of window sizes.")
		
	#Print an error if number of files incorrect.
	else:
		print("Error: number of files in input directory not equal to number of window sizes.")
	
"""
Match each input with the given window size to the nearest cluster for that window size.
Print out the region with its corresponding cluster to a BED file.
"""
def match_clusters(in_file_name, win, cluster_file_name, out_dir):
	
	#Open input and cluster files.
	in_file = open(in_file_name, "r")
	cluster_file = open(cluster_file_name, "r")
	
	#Open output file.
	out_file = open(out_dir + "map" + str(win) + ".bed", "w")
	
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
			[match, ambig] = match_datum(input, clusters)
			
			#Print match to BED file. Format is:
			#chrom	start	end	cluster_num	1 - ambiguity
			#Score is the opposite of the ambiguity metric.
			out_file.write("chr" + labels[0] + "\t" + labels[1] + "\t" + labels[2] + "\t" + str(match) + "\t" + str(1 - ambig) + "\n")

			#Read the next line in the file.			
			next_line = in_file.readline()
			
		#Print a message to the user.
		print("Files done for window size " + str(win))
	else:
		print("Only one cluster for window size " + str(win) + ". No annotation performed.")
	out_file.close()
	in_file.close()

"""
Find the closest match for the datum in the list of clusters.
Return an ambiguity metric which measures how close the datum
is to its nearest cluster as opposed to other clusters.
"""
def match_datum(datum, clusters):

	#Create array to hold distances between datum and clusters.
	max_crosscorr = 0
	match = 0
	crosscorr_list = []
	
	#For each cluster, determine the datum's distance from it.
	for i in range(len(clusters)):
		cluster = clusters[i]
		crosscorr_list.append(get_max_crosscorr(datum, cluster))
		if crosscorr_list[i] > max_crosscorr:
			max_crosscorr = crosscorr_list[i]
			match = i
	
	#Calculate the ambiguity metric
	ambig = get_ambiguity(crosscorr_list, datum)
			
	#The ambiguity metric is the ratio of the closest cluster's distance
	#to the distance of the second closest cluster.
	return [match, ambig]

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
def get_max_crosscorr(datum, cluster):

	#Since the datum is twice the size of the cluster, search along the datum and match at the best point.
	maximum = 0
	crosscorr = 0
	for delay in range(0, len(cluster)):
		try:
			crosscorr = co.get_crosscorr(datum, cluster, delay, 0.25, 0, False, False)
				
			#If the distance is the smallest so far, update.
			if crosscorr > maximum:
				maximum = crosscorr
		except ValueError:
			pass
	#Return the smallest distance
	return max(crosscorr, 0.01)
	
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