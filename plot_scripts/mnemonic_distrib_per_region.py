"""
This script finds the distribution of mnemonic
percentages across regions, for the following
mnemonic categories: Promoter, Enhancer, Weak, Repressor.
The following arguments are required:
1. The ChromHMM file containing the mnemonics
2. The region size
"""
import sys
import pandas as pd

def main():

    # Arguments
    chromhmm = pd.read_csv(open(sys.argv[1], "r"), sep = "\t", header = None)
    region_size = int(sys.argv[2])
    out = sys.argv[3]
    
    # Find percentages for each mnemonic type.
    promoter = {"1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"}
    enhancer = {"6_EnhG", "7_Enh", "12_EnhBiv"}
    repressor = {"13_ReprPC", "14_ReprPCWk"}
    weak = {"9_Het", "15_Quies"}
    percents = find_percents(chromhmm, region_size, promoter, enhancer, repressor, weak)
    
    # Save percentages.
    percents.to_csv(open(out, "w"), index = False)
    
"""
Build a data frame of percentages for each region.
"""
def find_percents(chromhmm, region_size, promoter, enhancer, repressor, weak):

    # Percentages
    promoter_percent = []
    enhancer_percent = []
    repressor_percent = []
    weak_percent = []
    starts = []
    ends = []
    chroms = []
    
    # Starting and ending chromosome positions.
    i = 0
    chrom = chromhmm.iloc[i,0]
    start = chromhmm.iloc[i,1]
    end = start + region_size
    
    # Get percentage for each region. Loop through as long as we haven't reached the
    # last position. Ignore sex chromosomes.
    while chrom != "chrX" and chrom != "chrY" and chrom != "chrM" and i < chromhmm.shape[0]:
    
        # Get percentages and append to list.
        if i % 10000 == 0:
            print("Percentages complete up to " + chrom + " " + str(start))
        [p, e, r, w, i] = get_percentages(chrom, chromhmm, start, end, promoter, enhancer, repressor, weak, i)
        if p != 0 or e != 0 or r != 0 or w != 0:
            promoter_percent.append(p)
            enhancer_percent.append(e)
            repressor_percent.append(r)
            weak_percent.append(w)
            starts.append(start)
            ends.append(end)
            chroms.append(chrom)
        
        # Get new region and chromosome values.
        if i < chromhmm.shape[0] and chromhmm.iloc[i,0] != chrom:
            chrom = chromhmm.iloc[i,0]
            start = 0
            end = region_size
        else:
            start = end
            end = end + region_size 
        
    # Return a data frame containing all percentages.
    percent_df = pd.DataFrame({"Chromosome": chroms, "Start": starts, "End": ends, "Promoter": promoter_percent, "Enhancer": enhancer_percent, "Repressor": repressor_percent, "Weak": weak_percent}, columns = ["Chromosome", "Start", "End", "Promoter", "Enhancer", "Repressor", "Weak"])
    return percent_df
    
"""
Build a data frame of percentages for each region.
"""
def get_percentages(chrom, chromhmm, start, end, promoter, enhancer, repressor, weak, i):
    percent_p = percent_e = percent_r = percent_w = 0
    
    # For each mnemonic, add to percentages.
    stop = False
    while i < chromhmm.shape[0] and chromhmm.iloc[i,0] == chrom and chromhmm.iloc[i,1] < end and not stop:
    
        # Get the size of the mnemonic (within the region).
        mnemonic_size = chromhmm.iloc[i,2] - chromhmm.iloc[i,1]
        if chromhmm.iloc[i,2] >= end:
            mnemonic_size = mnemonic_size - (chromhmm.iloc[i,2] - end)
        if chromhmm.iloc[i,1] < start:
            mnemonic_size = mnemonic_size - (start - chromhmm.iloc[i,1])
            
        # Add the mnemonic to the percentages.
        if chromhmm.iloc[i,3] in promoter:
            percent_p = percent_p + (mnemonic_size / (end - start))
        elif chromhmm.iloc[i,3] in enhancer:
            percent_e = percent_e + (mnemonic_size / (end - start))
        elif chromhmm.iloc[i,3] in repressor:
            percent_r = percent_r + (mnemonic_size / (end - start))
        elif chromhmm.iloc[i,3] in weak:
            percent_w = percent_w + (mnemonic_size / (end - start))

        # If we haven't reached the end of the region
        # yet, continue. If the end of this mnemonic runs
        # past the end of the region, stop.
        if chromhmm.iloc[i,2] < end:
            i = i + 1
        else:
            stop = True
            
    return [percent_p, percent_e, percent_r, percent_w, i]
    
if __name__ == "__main__":
    main()