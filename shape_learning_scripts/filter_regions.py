"""
Filter the training region file so that it excludes regions in promoters.
Required arguments:
1. The input region file
2. The BED file containing the promoter regions.
3. The output region file
"""
import sys
import os
import pandas as pd
import pickle as pkl
from tqdm import tqdm
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs
import numpy as np
import random

def main():

    #Read in the bed file and shape data.
    regions = pkl.load(open(sys.argv[1], "rb"))
    promoters = pd.read_csv(open(sys.argv[2], "rb"), sep = "\t", header = None)
    out_regions = sys.argv[3]
    
    # Sort the regions by start position.
    starts = []
    print("Extracting region start positions")
    for region in tqdm(regions):
        starts.append(region.start)
    print("Sorting regions by position")
    order = np.argsort(np.asarray(starts))
    regions = [regions[i] for i in order]
    
    # Loop through the annotations and delete anything that overlaps.
    i = 0
    for j in tqdm(range(promoters.shape[0] - 1)):
    
        # Until we have reached a region that ends when/after the promoter
        # starts, keep moving forward.
        while i < len(regions) and regions[i].end < promoters.iloc[j,1]:
            i = i + 1
        while i < len(regions) and ((regions[i].start >= promoters.iloc[j,1] and regions[i].start < promoters.iloc[j,2]) or (regions[i].end > promoters.iloc[j,1] and regions[i].end <= promoters.iloc[j,2])):
            regions.pop(i)
      
    # Shuffle and output the regions.
    random.shuffle(regions)
    pkl.dump(regions, open(out_regions, "wb"))
    
if __name__ == "__main__":
    main()