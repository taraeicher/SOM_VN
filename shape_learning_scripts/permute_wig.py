import numpy as np
import pandas as pd
import sys
import os
from numpy import random

"""
Permute the WIG signals.
"""
def main():
    #Read the WIG file, skipping the header lines, and extract the column with the signal intensity.
    #Permute the values in this column and create a new data file.
    #Write the new data to the new file.
    in_wig = pd.read_csv(sys.argv[1], delimiter = "\t", skiprows = [0,1])
    intensities = in_wig.iloc[:,1]
    int_perm = pd.DataFrame(np.random.permutation(intensities))
    new_data = pd.concat([in_wig.iloc[:,0], int_perm], axis=1)
    new_data.to_csv(sys.argv[2], sep='\t', header = False, index = False)
    
    #Write the headers.
    in_wig_2 = open(sys.argv[1], "r")
    header_1 = in_wig_2.readline()
    header_2 = in_wig_2.readline()
    in_wig_2.close()
    new_file = open(sys.argv[2], "r")
    new_file_data = new_file.read()
    new_file.close()
    actual_new_file = open(sys.argv[2] + "_tmp", "w")
    actual_new_file.write(header_1)
    actual_new_file.write(header_2)
    actual_new_file.write(new_file_data)
    os.rename(sys.argv[2] + "_tmp", sys.argv[2])

if __name__ == "__main__":
    main()