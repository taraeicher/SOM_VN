"""
Save a pickled region file as a CSV with a
column for each signal.
Required arguments:
1. The pickled region file
2. The location where you would like to save the CSV file.
"""
import sys
import pickle as pkl
import os
sys.path.append(os.path.abspath("../common_scripts"))
import region_defs

def main():

    #Read in the bed file and shape data.
    pickled_regions = pkl.load(open(sys.argv[1], "rb"))
    output_file = open(sys.argv[2], "w")
    
    # Loop through the regions and write them to a file.
    for region in pickled_regions:
        output_file.write(region.chromosome + "," + str(region.start) + "," + str(region.end) + "," + ",".join([str(signal) for signal in region.signals]) + "\n")
        
    # Close the file.
    output_file.close()
    
if __name__ == "__main__":
    main()