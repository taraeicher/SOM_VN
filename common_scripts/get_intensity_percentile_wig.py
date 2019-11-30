import wig_and_signal_utils as wsu
import numpy as np
import math
import sys
def main():
    # Input
    wig_file = open(sys.argv[1], "r")
    resolution = int(sys.argv[2])
    region_size = float(sys.argv[3])
    percentile = float(sys.argv[4])
    percentile_cutoff_file = open(sys.argv[5], "w")
    
    # Percentile
    percentile_cutoff = wsu.get_intensity_percentile(percentile, wig_file, resolution)
    percentile_cutoff_file.write(str(percentile_cutoff))
    percentile_cutoff_file.close()
    
if __name__ == "__main__":
    main()