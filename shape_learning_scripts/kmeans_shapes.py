"""
This script takes in a set of shapes and clusters them using K-means, where
K is determined using the Gap Statistic. This script depends on the gap-stat
package, which you can download from https://github.com/milesgranger/gap_statistic
using pip.

The script takes two arguments:
1. The file containing the shapes to be clustered.
2. The file where the clustered shapes should be stored.
"""
"""
https://github.com/minddrummer/gap/blob/master/gap/gap.py
from urllib import request

def download(url):
    filename = url.split('/')[-1]
    print("Downloading" + str(filename))
    with request.urlopen('http://python.org/') as f:
        data = f.read()
        f.close()
    with open(filename, 'wb') as myfile:
        myfile.write(data)

# get repository
download('https://github.com/minddrummer/gap/blob/master/gap/gap.py')
gap = imp.load_source('gap', '/users/PAS0272/osu5316/miniconda3/lib/python3.5/site-packages/gap/gap.py')
#gap_src.gap()

from setuptools import setup

setup(name='kmeans_shapes.py',
      version='0.1',
      author='Tara Eicher',
      dependency_links=['https://github.com/minddrummer/gap/blob/master/gap/gap.py'])
      
#from gap import gap
"""
from sklearn.cluster import KMeans
import numpy as np
import sys
import os
import imp
import pickle as pkl
import region_defs
import gap-stat

"""
Cluster the SOM using hierarchical clustering.
"""
def main():

    file_path = sys.argv[1]
    output_path = sys.argv[2]

    #Obtain the SOM shapes and cluster them.
    
    #Open the file containing the SOM shapes.
    #Add all shapes to the list of data to cluster.
    som_shapes = []
    if os.path.exists(file_path):
        shapes = pkl.load(open(file_path, 'rb'))
        for shape in shapes:
            som_shapes.append(shape.signals)
        
        #Compute the gap statistic and use it to find the best k-value.
        gaps, s_k, K = gap.gap_statistic(np.array(som_shapes), refs=None, B=10, K=range(1,len(som_shapes)), N_init = 10)
        best_k_value = gap.find_optimal_k(gaps, s_k, K)
        
        #Print message to user.
        print("Optimal K is " + str(best_k_value))
        
        #Perform k-means clustering.
        kmeans = KMeans(n_clusters=best_k_value, random_state=0).fit(np.array(som_shapes))
                        
        #Print cluster shapes from hierarchical clustering.
        file = open(output_path, 'w')
        print_shapes(kmeans.cluster_centers_, file)
            
    #Print message to user.
    print("Clustering complete.")
    
#Print out the cluster shapes
def print_shapes(cluster_centers, file):

    #Find the number of clusters.
    num_clusters = len(cluster_centers)
    
    #Print each cluster center.
    shapes = []
    for c in range(0, num_clusters):
        #Print the cluster center.
        shapes.append(region_defs.Shape(c, len(c), cluster_centers[c][val]))
    pkl.dump(shapes, file)   

if __name__ == "__main__":
    main()