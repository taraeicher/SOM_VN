"""
Given a CSV file with shape signals, save the shapes to a pickle file.
Required arguments:
1. The input CSV file
2. The output pickle file
"""
import sys
import pickle as pkl
import pandas as pd
import numpy as np
import os
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs

def main():

    #Read in the bed file and shape data.
    csv_shapes = pd.read_csv(open(sys.argv[1], "r"), header = None)
    output_file = sys.argv[2]
    
    # Loop through the shapes and build the pickled list.
    # Map count is not tracked in the CSV file, and it is not used
    # downstream, so we set it to 0.
    shapes = []
    for i in range(csv_shapes.shape[0]):
        shapes.append(region_defs.Shape(str(i), 0, np.asarray(csv_shapes.iloc[i])))
        print(shapes[i].name)
        
    # Save the pickled list.
    pkl.dump(shapes, open(output_file, "wb"))
    
if __name__ == "__main__":
    main()