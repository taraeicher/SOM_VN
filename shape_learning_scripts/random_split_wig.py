"""
Randomly split the WIG bins 50/50.
Required arguments:
1. The input region file.
2. The output region file for training.
3. The output region file for testing.
"""
import sys
import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import math

def main():

    #Read in the regions and permute.
    bins = pd.read_csv(open(sys.argv[1], "rb"), sep = "\t")
    bins = bins.iloc[np.random.permutation(bins.index)]
    bins.iloc[list(range(math.ceil(bins.shape[0] / 2)-1))].to_csv(open(sys.argv[2], "w"), sep = "\t", index = False, header = False)
    bins.iloc[list(range(math.ceil(bins.shape[0] / 2), bins.shape[0]-1))].to_csv(open(sys.argv[3], "w"), sep = "\t", index = False, header = False)
    
if __name__ == "__main__":
    main()