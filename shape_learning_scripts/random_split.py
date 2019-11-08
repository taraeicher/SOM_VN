"""
Randomly split the training regions 50/50.
Required arguments:
1. The input region file.
2. The output region file for training.
3. The output region file for testing.
"""
import sys
import os
import pickle as pkl
from tqdm import tqdm
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs
import random
import math

def main():

    #Read in the regions and permute.
    regions = pkl.load(open(sys.argv[1], "rb"))
    random.shuffle(regions)
    pkl.dump(regions[0:int(math.ceil((len(regions) / 2)))-1], open(sys.argv[2], "wb"))
    pkl.dump(regions[int(math.ceil((len(regions) / 2))):len(regions)-1], open(sys.argv[3], "wb"))
    
if __name__ == "__main__":
    main()