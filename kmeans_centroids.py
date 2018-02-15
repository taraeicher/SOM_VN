"""
Requires gap statistic. Need to import from https://github.com/minddrummer/gap/blob/master/gap/gap.py
"""
from sklearn.cluster import KMeans
import numpy as np
import sys
from gap import gap

"""
Cluster the SOM using hierarchical clustering.
"""
def main():

	num_windows = int(sys.argv[1])
	file_path = sys.argv[2]
	output_path = sys.argv[3]

	#For each file, obtain the SOM centroids and cluster them.
	for win in range(0, num_windows):
	
		#Open the file containing the SOM centroids.
		#Add all centroids to the list of data to cluster.
		som_centroids = []
		file = open(file_path + str(win), 'r')
		next_line = file.readline()
		while next_line:
			som_centroids.append([float(i) for i in next_line.split(",")])
			next_line = file.readline()
		
		#Compute the gap statistic and use it to find the best k-value.
		gaps, s_k, K = gap.gap_statistic(np.array(som_centroids), refs=None, B=10, K=range(1,len(som_centroids)), N_init = 10)
		bestKValue = gap.find_optimal_k(gaps, s_k, K)
		
		#Print message to user.
		print("Optimal K for window " + str(win) + " is " + str(bestKValue))
		
		#Perform k-means clustering.
		kmeans = KMeans(n_clusters=bestKValue, random_state=0).fit(np.array(som_centroids))
						
		#Print cluster centroids from hierarchical clustering.
		file = open(output_path + str(win), 'w')
		print_centroids(kmeans.cluster_centers_, som_centroids, file)
		
	#Print message to user.
	print("Clustering complete for all windows.")
	
#Print out the cluster centroids
def print_centroids(cluster_centers, som_centroids, file):

	#Find the number of clusters.
	num_clusters = len(cluster_centers)
	
	#For each cluster, find its mean.
	for c in range(0, num_clusters):
		#Print the cluster center.
		for val in range(len(cluster_centers[0])):
			file.write(str(cluster_centers[c][val]))
			if val < len(cluster_centers[0]) - 1:
				file.write(",")
		file.write("\n")

if __name__ == "__main__":
	main()