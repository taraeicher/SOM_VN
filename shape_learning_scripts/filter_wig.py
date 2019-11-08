"""
Filter the WIG file so that it excludes regions in promoters.
Required arguments:
1. The input WIG file
2. The BED file containing the promoter regions.
3. The output WIG file
"""
import sys
import os
import pandas as pd
from tqdm import tqdm
import numpy as np

def main():

    #Read in the bed file and shape data.
    wig = pd.read_csv(open(sys.argv[1], "rb"), sep = "\t", header = None)
    promoters = pd.read_csv(open(sys.argv[2], "rb"), sep = "\t", header = None, skiprows = 1)
    out_wig = sys.argv[3]
    chrom = sys.argv[4]

    # Loop through the regions and write them to a file.
    bins_overlapping_promoter = []
    for i in tqdm(np.where(promoters.iloc[:,0] == "chr" + chrom)[0][1:100]):
        where_starts_after_beginning = np.where(np.asarray(wig.iloc[:,1]) > promoters.iloc[i,1])[0]
        where_starts_before_end = np.where(np.asarray(wig.iloc[:,1]) <= promoters.iloc[i,2])[0]
        where_ends_after_beginning = np.where(np.asarray(wig.iloc[:,2]) >= promoters.iloc[i,1])[0]
        where_ends_before_end = np.where(np.asarray(wig.iloc[:,2]) < promoters.iloc[i,2])[0]
        where_starts_in_region = set(where_starts_after_beginning).intersection(set(where_starts_before_end))
        where_ends_in_region = set(where_ends_after_beginning).intersection(set(where_ends_before_end))
        bins_overlapping_promoter.extend(list(where_starts_in_region.union(where_ends_in_region)))
    wig = wig.drop(wig.index[bins_overlapping_promoter])
        
    # Close the file.
    wig.to_csv(open(out_wig, "w"), sep = "\t", index = False, header = False)
    
if __name__ == "__main__":
    main()