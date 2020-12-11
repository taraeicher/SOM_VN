#Import numpy for calculations
import numpy as np

#Import glob for obtaining the files.
import glob, os

#Import sys for obtaining command line args.
import sys
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
Match each input to its closest shape. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    input = sys.argv[1]
    shape = sys.argv[2]
    output_dir = sys.argv[3]
    wig_path = sys.argv[4]
    cutoff = float(sys.argv[5])

    #Open all input files.
    in_file = open(input, 'r')
    magnitude_file = open(shape, 'r')

    #Notify the user that files were found
    print("Using the input file:")
    print(in_file.name)

    #Notify the user that files were found
    print("Using the shape file:")
    print(magnitude_file.name)

    #Match the inputs to magnitudes and close the files.
    match_magnitudes_cutoff(in_file.name, magnitude_file.name, output_dir, wig_path, cutoff)
    in_file.close()
    magnitude_file.close()

"""
Match each input to the nearest shape.
Print out the region with its corresponding shape to a BED file.
"""
def match_magnitudes_cutoff(in_file_name, magnitude_file_name, out_dir, wig_name, cutoff):
    
    #Get threshold to use in printing.
    wig = open(wig_name, 'r')
    intensity = wsu.get_intensity_percentile(0.995, wig, 0)
    if intensity == 0:
        intensity = 0.1
    print("\n")
    wig.close()
    scale = 5.0 / intensity
    
    #Open input and shape files.
    in_file = open(in_file_name, "r")
    magnitude_file = open(magnitude_file_name, "r")
    
    #Open output files.
    out_file = open(out_dir, "w")
    out_clust = open(out_dir + "clust", "w")
    
    #Read in shape data.
    magnitudes = []
    magnitude_anno = []
    magnitude_name = []
    next_magnitude = magnitude_file.readline()
    counter = 0
    while next_magnitude:
        split_tabs = next_magnitude.split("\t")
        magnitude_anno.append(split_tabs[1])
        magnitude_name.append(split_tabs[0])
        magnitudes.append(int(split_tabs[2]))
        next_magnitude = magnitude_file.readline()
    magnitude_file.close()

    #Read in each line in the file and map it.
    if(len(magnitudes) > 1):
        next_line = in_file.readline()
        while(next_line):
        
            #Get the chromosome and position info to print to the output file.
            #Get the input signal to match with the magnitudes.
            inputStr = []
            labels = []
            split_line = next_line.split(",")
            labels = list(split_line[0:3])
            inputStr = list(split_line[3:len(split_line)])
            input = [float(i) for i in inputStr] 
            input_scaled = input * np.tile(scale, len(input))
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [match, ambig, out_str] = match_region(input_scaled, magnitudes, magnitude_anno)
            anno_label = magnitude_anno[match]
            anno_name = magnitude_name[match]

            #Print match to BED file. Format is:
            #chrom  start   end shape_num 1 - ambiguity
            #Score is the opposite of the ambiguity metric.
            out_file.write("chr" + labels[0] + "\t" + labels[1] + "\t" + labels[2] + "\t" + anno_label + "\t" + str(1 - ambig) + "\t" + out_str + "\n")
            out_clust.write("chr" + labels[0] + "\t" + labels[1] + "\t" + labels[2] + "\t" + anno_name + "\t" + str(1 - ambig) + "\t" + out_str + "\n")

            #Read the next line in the file.            
            next_line = in_file.readline()
            
        #Print a message to the user.
        print("Files done")
    else:
        print("Only one magnitude. No annotation performed.")
    out_file.close()
    out_clust.close()
    in_file.close()

"""
Find the closest match for the input in the list of magnitudes.
"""
def match_region(region, magnitudes, annotations):

    #Create array to hold distances between region and magnitudes.
    match = 0
    dist_list = []
    max_region = np.max(region)
    min_dist = np.max(magnitudes)
    
    #For each shape, determine the region's distance from it.
    for i in range(len(magnitudes)):
        magnitude = magnitudes[i]
        if annotations[i] != "Unknown":
            dist = abs(max_region - magnitude)
            dist_list.append(dist)
            if dist_list[len(dist_list) - 1] < min_dist:
                min_dist = dist_list[len(dist_list) - 1]
                match = i

    ambig = get_ambiguity(dist_list)
    
    #Write the magnitude values to a separate file.
    #This file will be used for shape validity analysis.
    comma = ","
    out_str = comma.join(str(e) for e in region)
            
    #The ambiguity metric is the ratio of the closest shape's distance
    #to the distance of the second closest shape.
    return [match, ambig, out_str]
    
"""
Compute the ambiguity metric. It is equivalent to the smallest distance
from the region to a magnitude, divided by the second smallest
distance from the region to a magnitude.
""" 
def get_ambiguity(l):
    list_array = np.asarray(l)
    min_idx = np.argmin(list_array)
    list_min_removed = np.delete(list_array, min_idx)
    min_idx_2 = np.argmin(list_min_removed)
    ambiguity = list_array[min_idx] / list_min_removed[min_idx_2]
    return ambiguity

if __name__ == "__main__":
    main()