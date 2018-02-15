from scipy import cluster
import numpy as np
import sys
import math

"""
Cluster the SOM using hierarchical clustering.
"""
def main():

	num_windows = int(sys.argv[1])
	file_path = sys.argv[2]
	output_path = sys.argv[3]

	#For each file, obtain the SOM centroids and cluster them.
	for win in range(0, num_windows):
	
		#Print message to user.
		print("Computing hierarchical clustering for window " + str(win))
	
		#Open the file containing the SOM centroids.
		#Add all centroids to the list of data to cluster.
		som_centroids = []
		file = open(file_path + str(win), 'r')
		next_line = file.readline()
		while next_line:
			som_centroids.append([float(i) for i in next_line.split(",")])
			next_line = file.readline()
		
		#Get the linkages for the hierarchical clusters.
		linkage_matrix = cluster.hierarchy.linkage(som_centroids)
		
		#Use the RMSSTD to determine which inconsistency metric cutoff is optimal.
		#We want to find the location where there is an "elbow" (i.e. a significant
		#change in value) and choose the cutoff where the elbow occurs.
		elbow_slope = 0
		prev_rmsstd = 0
		prev_t = 0
		opt_t = 0
		for t_int in range(10, 30):
			t = float(t_int) / 10
			cluster_list = cluster.hierarchy.fcluster(linkage_matrix, t, depth = 10)
			rmsstd = get_rmsstd(cluster_list, som_centroids)
			if t > 1 and rmsstd - prev_rmsstd > elbow_slope:
				elbow_slope = rmsstd - prev_rmsstd
				opt_t = prev_t
			prev_rmsstd = rmsstd
			prev_t = t
		opt_cluster_list = cluster.hierarchy.fcluster(linkage_matrix, opt_t, depth = 10)
		print(opt_t)
						
		#Print cluster centroids from hierarchical clustering.
		file = open(output_path + str(win), 'w')
		print_centroids(opt_cluster_list, som_centroids, file)
		
	#Print message to user.
	print("Clustering complete for all windows.")
	
#Calculate the root mean square standard deviation for the set of clusters.
# A low RMSSTD indicates a good clustering, so we want to minimize this.
def get_rmsstd(cluster_list, som_centroids):
	num_vars = len(som_centroids[1])
	numerator = 0
	denominator = 0
	
	#Find the number of clusters.
	num_clusters = max(cluster_list)
	
	#For each cluster, find its count and sum of squares.
	for c in range(1, num_clusters + 1):
		cluster = []
		
		#Find each element of the cluster.
		for i in range(len(cluster_list)):
			if cluster_list[i] == c:
				cluster.append(som_centroids[i])
		
		#Add the calculations for this cluster to the total.
		clust_array = np.array(cluster)
		mean = np.mean(clust_array, 0)
		denominator += num_vars * len(cluster) - 1
		for centroid in cluster:
			numerator += np.sum(np.power(np.subtract(centroid, mean), 2))

	#Return the square root.
	return math.sqrt(float(numerator) / denominator)
	
#Print out the cluster centroids
def print_centroids(cluster_list, som_centroids, file):

	#Find the number of clusters.
	num_clusters = max(cluster_list)
	
	#For each cluster, find its mean.
	for c in range(1, num_clusters + 1):
		cluster = []
		#Find each element of the cluster.
		for i in range(len(cluster_list)):
			if cluster_list[i] == c:
				cluster.append(som_centroids[i])
				
		#Find the mean.
		clust_array = np.array(cluster)
		mean = np.mean(clust_array, 0)
		
		#Print the mean.
		for val in range(len(mean)):
			file.write(str(mean[val]))
			if val < len(mean) - 1:
				file.write(",")
		file.write("\n")

if __name__ == "__main__":
	main()