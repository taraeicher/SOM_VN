import numpy as np
import pandas as pd
import sys
import os
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu
BIN_SIZE = 50
"""
Given the WIG file, predict using RPKM percentile and write to a BED file.
We first annotate each bin, then merge bins with the same annotation.
"""
def main():

    #Get input, output path, and chromosome.
    wig = sys.argv[1]
    bin_out = sys.argv[2]
    out_merged = sys.argv[3]
    chrom = sys.argv[4]
    
    #First pass: Annotate each bin and output the file.
    annotate_all_bin_bed(wig, bin_out, chrom)
    merge_bed(bin_out, out_merged, chrom)

"""
Annotate each 50 bp bin individually using RPKM cutoffs.
"""
def annotate_all_bin_bed(wig_f, out_f, chrom):

    #Open the input and output.
    wig = open(wig_f, "r")
    out = open(out_f, "w")
    
    #We choose the 97th and the 90th percentiles for promoter and enhancer, respectively.
    promoter_cutoff = wsu.get_intensity_percentile(0.97, wig, 0)
    enhancer_cutoff = wsu.get_intensity_percentile(0.90, wig, 0)

    #Loop through the WIG file and annotate each region. Skip the first two lines, which are the header.
    junk = wig.readline()
    junk = wig.readline()
    bin = wig.readline().split("\t")
    while len(bin) > 1:
        
        #Output the starting and ending positions.
        end = str(int(bin[0]) + BIN_SIZE)
        out.write("chr" + chrom + "\t" + bin[0] + "\t" + end + "\t")
        
        #Get the annotation.
        rpkm = float(bin[1])
        if rpkm >= promoter_cutoff:
            out.write("Promoter\n")
        elif rpkm >= enhancer_cutoff:
            out.write("Enhancer\n")
        else:
            out.write("Weak\n")
        
        #Get the next line.
        bin = wig.readline().split("\t")
            
    #Close the input and output.
    wig.close()
    out.close()
    
"""
Given a BED file, merge all adjacent regions with the same annotation.
"""
def merge_bed(bed_f, out_f, chrom):

    #Open the input and output files.
    bed = open(bed_f, "r")
    out = open(out_f, "w")
    
    #Loop through all lines in the BED file, merging them whenever the annotations match.
    bin = bed.readline().split("\t")
    first_line = True
    running_bin = bin
    while len(bin) > 1:
        
        #If the next line has the same annotation as the current line and immediately
        #follows it, add it to the running bin.

        if bin[3].rstrip() == running_bin[3].rstrip() and int(bin[1]) == int(running_bin[2]):
            running_bin[2] = bin[2]
            
        #Otherwise, output the current running bin and start a new one.
        elif not first_line:
            out.write("\t".join(running_bin))
            running_bin = bin
            
        #Get the next bin.
        bin = bed.readline().split("\t")
        first_line = False
        
    #Close the input and output.
    bed.close()
    out.close()

if __name__ == "__main__":
    main()