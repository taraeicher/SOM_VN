import pandas as pd
import numpy as np
import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#Code from https://stackoverflow.com/questions/48036593/is-my-python-implementation-of-the-davies-bouldin-index-correct?rq=1
def DaviesBouldin(X, centroids, labels):
    n_cluster = len(np.bincount(labels))
    cluster_k = [X[labels == k] for k in range(n_cluster)]
    centroids = [np.mean(k, axis = 0) for k in cluster_k]
    variances = [np.mean([np.linalg.norm(p - centroids[i]) for p in k]) for i, k in enumerate(cluster_k)]
    db = []
    nans = set()

    for i in range(n_cluster):
        for j in range(n_cluster):
            #Add the Davies-Bouldin metric for this pair to the list.
            if j != i and not math.isnan(variances[i]) and not math.isnan(variances[j]):
                db.append((variances[i] + variances[j]) / np.linalg.norm(centroids[i] - centroids[j]))
            #If variance is nan, then there were no elements in the cluster. Remove it.
            elif math.isnan(variances[i]):
                nans.add(i)
            elif math.isnan(variances[j]):
                nans.add(j)
                
    #Subtract the number of clusters with no elements from the total number of clusters.
    n_cluster = n_cluster - len(nans)

    return(np.max(db) / n_cluster)

def main():
    bed_dir = sys.argv[1]
    centroids_dir = sys.argv[2]
    clusters_dir = sys.argv[3]
    heat_path = sys.argv[4]
    
    #Compute the index for each window size and each chromosome.
    windows = [2, 3, 4, 5, 6]
    window_sizes = ['2000', '4000', '8000', '16000', '32000']
    chromosomes = ['4', '8', '9', '13', '14', '15', '16', '17', '18', '19', '20', '21']
    indices = np.ones((len(windows), len(chromosomes)))
    for i in range(0, len(windows)):
        for j in range(0, len(chromosomes)):
            bed = pd.read_csv(bed_dir + "anno" + chromosomes[j] + ".map" + str(windows[i]) + ".bed", sep = "\t", header = None)
            centroids = np.loadtxt(open(centroids_dir + "chrom" + chromosomes[j] + "som_centroid" + str(windows[i]), "rb"), delimiter=",")
            clusters = np.loadtxt(open(bed_dir + "clusters_anno" + chromosomes[j] + ".map" + str(windows[i]), "rb"), delimiter=",")
            labels = bed.loc[:,3]
            indices[i,j] = DaviesBouldin(clusters, centroids, labels)
    
    #Plot the heatmap.
    plot_heatmap(indices.transpose(), window_sizes, chromosomes, heat_path)
    
#Plot a heatmap of total counts and a stacked bar plot of annotation type.
def plot_heatmap(values, xnames, ynames, heat_path):
    
    #Plot a heatmap of total counts.
    f = plt.figure(1)
    ax = plt.axes()
    heatmap = sns.heatmap(values, cbar=True, cmap="binary_r", fmt="d", vmin = 0, vmax = 1, xticklabels = xnames, yticklabels = ynames, cbar_kws={'label': 'Davies-Bouldin Index Value'})
    ax.set(xlabel = "Window Size (bp)", ylabel = "Chromosome")
    plt.yticks(rotation = 0)
                       
    fig = heatmap.get_figure()
    fig.savefig(heat_path + ".png")

if __name__ == "__main__":
    main()