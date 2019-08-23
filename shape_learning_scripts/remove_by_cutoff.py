"""
This script takes in shapes learned by the SOM and removes those with
mapping counts below a cutoff. The filtered shapes are saved in a new file.
"""
import sys
import os
import pickle as pkl
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs
from tqdm import tqdm

"""
Remove all shapes that did not have at least *cutoff* regions mapping to them in the last iteration.
"""
def main():
    file_path = sys.argv[1] #file path for obtaining clusters and counts
    cutoff = int(sys.argv[2])   #count cutoff
    output_path = sys.argv[3]   #output clusters

    #Obtain the SOM shapes and cluster them.

    #Print message to user.
    print("Filtering grid by cutoff")

    #Open the file containing the SOM shapes.
    #Add all shapes to the list of data to cluster.
    filtered_shapes = []
    if os.path.exists(file_path):
        shapes = pkl.load(open(file_path, 'rb'))
        for shape in tqdm(shapes):
            if shape.mapped_count >= cutoff:
                filtered_shapes.append(shape)
                
        #Print all shapes.
        pkl.dump(filtered_shapes, open(output_path, 'wb'))
            
if __name__ == "__main__":
    main()
    