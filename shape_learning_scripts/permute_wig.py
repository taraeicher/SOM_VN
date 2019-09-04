"""
Permute the WIG signals. Note that we only permute on non-zero bins
and single bins with zeros, not on long stretches of zeros. The following
arguments are required:
1. The WIG directory.
2. The directory where the permuted WIG should be stored.
3. The bin size
"""
import numpy as np
import pandas as pd
import sys
import os
from numpy import random
from tqdm import tqdm

def main():
    # arguments
    in_wig_dir = sys.argv[1]
    out_wig_dir = sys.argv[2]
    
    # Create the output directory.
    if not os.path.exists(out_wig_dir):
        os.makedirs(out_wig_dir)
    
    for file in os.listdir(in_wig_dir):
        print("Permuting " + file)
        in_wig = pd.read_csv(in_wig_dir + "/" + file, delimiter = "\t", header = None)
        perm_df = permute(in_wig)
        perm_df.to_csv(out_wig_dir + "/" + file, sep='\t', header = False, index = False)
        
def permute(wig):

    # Permute all rows in the WIG file.
    permuted = pd.DataFrame(np.random.permutation(wig))

    # Get the sizes of all bins (may not be equal).
    diffs = permuted.iloc[:,2] - permuted.iloc[:,1]
    start = np.zeros(permuted.shape[0], dtype = "int")
    end = np.zeros(permuted.shape[0], dtype = "int")
    
    # Built the new start and end positions, retaining the
    # bin size.
    start_pos = 0
    for bin in tqdm(range(permuted.shape[0])):
        start[bin] = int(start_pos)
        end[bin] = int(start_pos + diffs.iloc[bin])
        start_pos = start_pos + diffs.iloc[bin]

    # Integrate the new positions into the data frame.
    start_df = pd.DataFrame(start)
    end_df = pd.DataFrame(end)
    permuted.iloc[:,1] = start_df.values
    permuted.iloc[:,2] = end_df.values
        
    return permuted
    
    
    
    
    
    
if __name__ == "__main__":
    main()