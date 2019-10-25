"""
This script splits a WIG file into regions with a given region size
and margin size. It takes in the following required parameters:
1. A WIG file containing RPKM signal for a single chromosome.
2. The resolution of the WIG file (base pairs per bin)
3. The chromosome (do not include 'chr').
4. The name of the file where you wish to output the regions.
5. The total region size (in bp).
6. The margin size (in bp)
7. The factor to use in computing the optimal shift.
8. The signal percentile to use as the cutoff for crossing count.
9. The name of the file where you wish to output the regions after shifting.
10. The file where you wish to store the sum of all median signals after shifting.
11. The file where you wish to store the percentile cutoff used.
"""
import sys
import region_defs
import math
import pickle as pkl
import random
import wig_and_signal_utils as wsu
from tqdm import tqdm

"""
main function
"""
def main():

    # Get arguments.
    wig_file = open(sys.argv[1], "r")
    resolution = int(sys.argv[2])
    chrom = sys.argv[3]
    output_file = open(sys.argv[4], "wb")
    region_size = float(sys.argv[5])
    margin = float(sys.argv[6])
    factor = float(sys.argv[7])
    percentile = float(sys.argv[8])
    output_file_shifted = open(sys.argv[9], "wb")
    crossing_counts_file = open(sys.argv[10], "w")
    percentile_cutoff_file = open(sys.argv[11], "w")
    
    # Build list of regions and save.
    regions = build_region_list(wig_file, resolution, region_size, chrom, margin)
    wig_file.close()
    random.shuffle(regions)
    
    # Shift the regions.
    wig_file = open(sys.argv[1], "r")
    [shifted_regions, crossings_counts, percentile_cutoff] = shift_all_regions(regions, int(region_size / resolution), percentile, factor, wig_file, resolution)
    wig_file.close()
    
    # Save all files.
    pkl.dump(regions, output_file)
    output_file.close()
    pkl.dump(shifted_regions, output_file_shifted)
    output_file_shifted.close()
    for crossing_count in crossings_counts:
        crossing_counts_file.write("%f\n" % crossing_count)
    crossing_counts_file.close()
    percentile_cutoff_file.write(str(percentile_cutoff))
    percentile_cutoff_file.close()

    # Notify user that process has completed.
    print("Input generated for chromosome " + chrom + " for margin size " + str(margin) + " and factor " + str(factor))

"""
For each region in the WIG file, shift it to its best
representation. Return the list of all shifted regions.
"""
def shift_all_regions(regions, region_size, percentile, factor, wig_file, bin_size):
    threshold = wsu.get_intensity_percentile(percentile, wig_file, bin_size)
    
    shifted_regions = []
    all_crossings = []
    print("Shifting regions:")
    for region in tqdm(regions):
       shifted_region = region_defs.Shifted_Region(region, region_size, factor, threshold)
       shifted_regions.append(shifted_region)
       all_crossings.append(shifted_region.crossings)
    return [shifted_regions, all_crossings, threshold]

"""
Build the list of all regions in the WIG file.
"""
def build_region_list(wig_file, resolution, region_size, chrom, margin):
    #Variables that will be used when traversing file.     
    interval_start = 0
    interval_end = 0
    chrom_pos = 0
    input = ""
    signal = 0.0
    region_size = region_size + 2 * margin
    region_list = []
    current_region = None

    #Extract the first signal.
    #Do not continue if the end of the file has been reached.
    input = wig_file.readline()
    signal = float(input.split("\t")[3])
    chrom_pos = float(input.split("\t")[1])
    chrom_pos_end = float(input.split("\t")[2])
    add_status = region_defs.Region.REGION_COMPLETED_WITH_SIGNAL
    while input is not None and input != "":

        # If the region was completed with the input signal, add it to the list
        # and start a new one.
        if add_status == region_defs.Region.REGION_COMPLETED_WITH_SIGNAL:

            # Add the region to the list.
            if current_region != None:
                region_list.append(current_region)
                #print(current_region.signals)
                
            # Get info needed for the region.
            interval_start = chrom_pos
            interval_end = interval_start + region_size
            
            # Create the region and add the signal.
            current_region = region_defs.Region(interval_start, interval_end, chrom, region_size / resolution)
            add_status = current_region.add_signal(chrom_pos, chrom_pos_end, signal, resolution)
            
            # Read the next line.
            input = wig_file.readline()
            if input is not None and input != "":
                signal = float(input.split("\t")[3])
                chrom_pos = float(input.split("\t")[1])
                chrom_pos_end = float(input.split("\t")[2])
            
        # If the region was completed with zeros, add it to the list and
        # start a new region closest to the current signal.
        elif add_status == region_defs.Region.REGION_COMPLETED_WITH_ZEROS:
        
            # Add the region to the list.
            if current_region != None:
                region_list.append(current_region)
                #print(current_region.signals)
                
            # Find the closest starting point to the current signal.
            interval_start = find_closest_starting_point(interval_start, region_size, chrom_pos)
            interval_end = interval_start + region_size
            current_region = region_defs.Region(interval_start, interval_end, chrom, region_size / resolution)
            add_status = current_region.add_signal(chrom_pos, chrom_pos_end, signal, resolution)
            
            # Read the next line.
            input = wig_file.readline()
            if input is not None and input != "":
                signal = float(input.split("\t")[3])
                chrom_pos = float(input.split("\t")[1])
                chrom_pos_end = float(input.split("\t")[2])
                
        # If the region was not completed, keep adding new signal.
        elif add_status == region_defs.Region.REGION_INCOMPLETE:
            add_status = current_region.add_signal(chrom_pos, chrom_pos_end, signal, resolution)
            
            # Read the next line.
            input = wig_file.readline()
            if input is not None and input != "":
                signal = float(input.split("\t")[3])
                chrom_pos = float(input.split("\t")[1])
                chrom_pos_end = float(input.split("\t")[2])
          
        # If the region was completed but there is still more signal in the bin,
        # create a new interval starting where we left off.
        elif add_status == region_defs.Region.REGION_SIGNAL_OVERFLOW:
        
            # Add the region to the list.
            if current_region != None:
                region_list.append(current_region)
                #print(current_region.signals)
                
            # Start a new region beginning with the current bin.
            interval_start = interval_end
            interval_end = interval_start + region_size
            current_region = region_defs.Region(interval_start, interval_end, chrom, region_size / resolution)
            add_status = current_region.add_signal(interval_start, chrom_pos_end, signal, resolution)
        
    return region_list
    
"""
Given the previous starting point, the next position, and the region
size, find a suitable starting point near the next position. The starting
point must be a multiple of the region size.
"""
def find_closest_starting_point(old_start, region_size, position):
    
    regions_away = math.floor((position - old_start) / region_size)
    new_start = old_start + regions_away * region_size
    return new_start

if __name__ == "__main__":
    main()