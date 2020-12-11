#Import sys for obtaining command line args.
import sys
import numpy as np
from scipy.stats import ranksums

"""
For each datum in each of the window size files, match it to its
closest cluster in the cluster files. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    #Window size list and list of chromosomes are variable.
    #Allow user to choose number of iterations, initial learning
    #rate, and initial neighborhood size.
    #20, 50, 100, 250,
    input_dir = sys.argv[1]
    input_dir_rand = sys.argv[2]

    #Get all cross-correlations.
    chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    cc_prom = []
    cc_enh = []
    cc_weak = []
    rand_cc_weak = []
    for c in chroms:

        #Get all cross-corr for the given annotation in that file.
        in_file = open(input_dir + c, 'r')
        in_file_rand = open(input_dir_rand + c, 'r')
        cc_prom.extend(get_crosscorr_anno(in_file, "Promoter"))
        cc_enh.extend(get_crosscorr_anno(in_file, "Enhancer"))
        cc_weak.extend(get_crosscorr_anno(in_file, "Weak"))
        rand_cc_weak.extend(get_crosscorr_anno(in_file_rand,"Weak"))
        in_file.close()
        in_file_rand.close()
   
    #Plot all cross-correlations.
    print("Promoter result")
    print(ranksums(cc_prom, rand_cc_weak))
    print("Enhancer result")
    print(ranksums(cc_enh, rand_cc_weak))
    print("Weak result")
    print(ranksums(cc_weak, rand_cc_weak))

"""
Match each input with the given window size to the nearest cluster for that window size.
Print out the region with its corresponding cluster to a BED file.
"""
def get_crosscorr_anno(in_file, annotation):

    #Read in each line in the file and add to the cross-correlation list if it is the correct annotation.
    cross_corrs = []
    next_line = in_file.readline()
    while(next_line):

        #Get the chromosome and position info to print to the output file.
        #Get the input signal to match with the clusters.
        split_line = next_line.split("\t")
        anno = split_line[3]
        cross_corr = float(split_line[4])

        #Match the data to the nearest cluster and obtain the match and the ambiguity metric.
        if anno == annotation:
            cross_corrs.append(cross_corr)

        #Read the next line in the file.
        next_line = in_file.readline()

    #Go bac to the beginning of the file and return the cross-corr list.
    in_file.seek(0)
    return cross_corrs

if __name__ == "__main__":
    main()
