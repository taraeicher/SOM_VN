import numpy as np
import pandas as pd
import sys
from numpy import random
"""
Permute the ChromHMM annotations. Annotated regions stay the same,
but the labels (e.g. Quiesc, Enh, etc) are rearranged.
"""
def main():
    #Read the input file and extract the column with the annotations.
    in_bed = pd.read_csv(sys.argv[1], delimiter = "\t", header = None, index_col = False)
    unique_chroms = in_bed.iloc[:,0].unique()
    for chrom in unique_chroms:
        chr_bed = in_bed.iloc[np.where(in_bed.iloc[:,0] == chrom)[0]].reindex()
        annotations = chr_bed.iloc[:,3]

        #Get starting position.
        start = int(chr_bed.iloc[0,1])
        
        #Compute length of all regions.
        lengths = np.asarray([int(i) for i in chr_bed.iloc[:,2]]) - np.asarray([int(i) for i in chr_bed.iloc[:,1]])
        
        #Permute the values in this column.
        indices = range(0, len(annotations))
        perm = np.random.permutation(indices)
        annotations_perm = np.asarray(annotations.iloc[perm])
        lengths_perm = lengths[perm]
        
        #Propagate lengths forward to get new indices.
        starts_perm = np.zeros(len(lengths_perm))
        ends_perm = np.zeros(len(lengths_perm))
        starts_perm[0] = start
        ends_perm[0] = start + lengths_perm[0]
        for i in range(1, len(starts_perm)):
            starts_perm[i] = ends_perm[i-1]
            ends_perm[i] = starts_perm[i] + lengths_perm[i]

        #Write the new data to the new file. Must do one line at a time, as pandas concat automatically sorts.
        new_data = pd.DataFrame({"Chrom": chr_bed.iloc[:,0], "Start": [str(int(i)) for i in starts_perm], "End": [str(int(i)) for i in ends_perm], "Annotation": annotations_perm}, columns = ["Chrom", "Start", "End", "Annotation"])
        print(new_data)
        new_data.to_csv(sys.argv[2] + "_" + chrom + ".bed", sep='\t', header = False, index = False)

if __name__ == "__main__":
    main()