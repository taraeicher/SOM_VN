#Import sys for obtaining command line args.
import sys
import os

"""
Remove all centroids that did not have at least *cutoff* regions mapping to them in the last iteration.
"""
def main():
    file_path = sys.argv[1] #file path for obtaining clusters and counts
    cutoff = int(sys.argv[2])   #count cutoff
    output_path = sys.argv[3]   #output clusters

    #Obtain the SOM centroids and cluster them.

    #Print message to user.
    print("Filtering grid by cutoff")

    #Open the file containing the SOM centroids.
    #Add all centroids to the list of data to cluster.
    som_centroids = []
    counts = []
    if os.path.exists(file_path):
        file = open(file_path, 'r')
        count_file = open(file_path + "_counts", 'r')
        next_line = file.readline()
        next_count = count_file.readline()
        while next_line and next_count:
            som_centroids.append(next_line)
            counts.append(float(next_count))
            next_line = file.readline()
            next_count = count_file.readline()
            
        #Remove all centroids that don't meet the criteria.
        j = 0
        while j < len(som_centroids):
            if counts[j] < cutoff:
                del som_centroids[j]
                del counts[j]
            else:
                j += 1
                
        #Print all centroids.
        out_file = open(output_path, 'w')
        for centroid in som_centroids:
            out_file.write(centroid)
            
if __name__ == "__main__":
    main()
    