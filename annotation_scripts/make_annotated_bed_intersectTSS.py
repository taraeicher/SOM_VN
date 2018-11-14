#Import numpy for calculations
import numpy as np
import pandas as pd

#Import glob for obtaining the files.
import glob, os

#Import sys for obtaining command line args.
import sys
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
For each datum in each of the window size files, match it to its
closest shape in the shape files. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    input = sys.argv[1]
    shape = sys.argv[2]
    output_dir = sys.argv[3]
    wig_path = sys.argv[4]
    cutoff = float(sys.argv[5])
    tss_file = sys.argv[6]

    #Create grids for labeling SOM nodes with shape indices.
    learned_shapes = []
    learned_shapes_counts = []
    learned_shapes.append(list())
    learned_shapes_counts.append(list())
    shape_names = [[]]
    shape_count = 0
    
    #Open all input files.
    in_file = open(input, 'r')
    shape_file = open(shape, 'r')

    #Notify the user that files were found
    print("Using the input file:")
    print(in_file.name)

    #Notify the user that files were found
    print("Using the shape file:")
    print(shape_file.name)
    
    #Load TSS file
    tss = pd.read_csv(tss_file, header = None, sep = " ")

    #Match the inputs to shapes and close the files.
    match_shapes_cutoff(in_file.name, shape_file.name, output_dir, wig_path, cutoff, tss)
    in_file.close()
    shape_file.close()

"""
Match each input with the given window size to the nearest shape for that window size.
Print out the region with its corresponding shape to a BED file.
"""
def match_shapes_cutoff(in_file_name, shape_file_name, out_dir, wig_name, cutoff, tss):
    
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
    shape_file = open(shape_file_name, "r")
    
    #Open output files.
    out_file = open(out_dir, "w")
    out_clust = open(out_dir + "clust", "w")
    
    #Read in shape data.
    shapes = []
    shape_anno = []
    shape_name = []
    next_shape = shape_file.readline()
    counter = 0
    while next_shape:
        split_tabs = next_shape.split("\t")
        shape_anno.append(split_tabs[1])
        shape_name.append(split_tabs[0])
        shapes.append([float(i) for i in split_tabs[2].split(",")])
        next_shape = shape_file.readline()
    shape_file.close()
    
    #Read in each line in the file and map it.
    if(len(shapes) > 1):
        next_line = in_file.readline()
        while(next_line):
        
            #Get the chromosome and position info to print to the output file.
            #Get the input signal to match with the shapes.
            inputStr = []
            labels = []
            split_line = next_line.split(",")
            labels = list(split_line[0:3])
            inputStr = list(split_line[3:len(split_line)])
            input = [float(i) for i in inputStr] 
            input_scaled = input * np.tile(scale, len(input))
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [match, ambig, crosscorr, out_str] = match_datum(input_scaled, shapes, shape_anno, in_tss(tss, labels))
            
            #If the cross-correlation is within the threshold, assign the label as the correct label for the region.
            #Otherwise, assign it as unknown.
            anno_label = "Unknown"
            anno_name = "Unknown"
            if crosscorr >= cutoff and match != -1:
                anno_label = shape_anno[match]
                anno_name = shape_name[match]

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
        print("Only one shape. No annotation performed.")
    out_file.close()
    out_clust.close()
    in_file.close()
    
"""
Check whether the region overlaps a transcription start site, where the
list of transcription start sites is given.
"""
def in_tss(tss_data, region_info):
    comes_before = np.where(tss_data[1] <= int(region_info[2]))[0]
    comes_after = np.where(tss_data[2] >= int(region_info[1]))[0]
    overlaps = set(comes_before).intersection(comes_after)
    return(len(overlaps) != 0)

"""
Find the closest match for the datum in the list of shapes.
Return an ambiguity metric which measures how close the datum
is to its nearest shape as opposed to other shapes.
"""
def match_datum(datum, shapes, annotations, keep_promoter):

    #Create array to hold distances between datum and shapes.
    max_crosscorr = 0
    match = 0
    opt_delay = 0
    crosscorr_list = []
    
    #For each shape, determine the datum's distance from it.
    for i in range(len(shapes)):
        shape = shapes[i]
        #Only consider promoters if we are inside a TSS.
        if annotations[i] != "Unknown" and (annotations[i] != "Promoter" or keep_promoter):
            crosscorr, d = get_max_crosscorr(datum, shape)
            crosscorr_list.append(crosscorr)
            if crosscorr_list[len(crosscorr_list) - 1] > max_crosscorr:
                max_crosscorr = crosscorr_list[len(crosscorr_list) - 1]
                match = i
                opt_delay = d

    #Calculate the ambiguity metric. Only do so if not all cross-correlations are zero.
    #If they are all zero, set ambiguity to 1 (completely ambiguous) and do not assign
    #it to a shape.
    ambig = 1.0
    if not max_crosscorr == 0:
        ambig = get_ambiguity(crosscorr_list, datum)
    else:
        match = -1
    
    #Write the shape values to a separate file.
    #This file will be used for shape validity analysis.
    
    comma = ","
    out_str = comma.join(str(e) for e in datum)
            
    #The ambiguity metric is the ratio of the closest shape's distance
    #to the distance of the second closest shape.
    return [match, ambig, max_crosscorr, out_str]

"""
Compute the ambiguity metric. It is equivalent to the smallest distance
from the datum to a shape centroid, divided by the second smallest
distance from the datum to a shape centroid.
""" 
def get_ambiguity(list, datum):
    list_array = np.asarray(list)
    max_idx = np.argmax(list_array)
    list_max_removed = np.delete(list_array, max_idx)
    max_idx_2 = np.argmax(list_max_removed)
    ambiguity = list_max_removed[max_idx_2] / list_array[max_idx]
    return ambiguity

"""
Find the minimum distance between the datum and shifted shape.
"""
def get_max_crosscorr(datum, shape):

    #Since the datum is twice the size of the shape, search along the datum and match at the best point.
    maximum = 0
    crosscorr = 0
    opt_delay = 0
    step = len(shape) / 10
    for i in range(0, 10):
        try:
            delay = i * step
            crosscorr = wsu.get_crosscorr(datum, shape, int(delay), 0.5, 0, False, False, 0)
  
            #If the distance is the smallest so far, update.
            if crosscorr > maximum:
                maximum = crosscorr
                opt_delay = delay
        except ValueError:
            pass
    #Return the smallest distance
    return [crosscorr, delay]
    
"""
Find the minimum distance between the datum and shifted shape.
"""
def get_min_shifted_distance(datum, shape):

    #Since the datum is twice the size of the shape, search along the datum and match at the best point.
    min_dist = sum(datum)
    for delay in range(0, len(shape)):
        datum_array = np.asarray(datum[delay:len(shape) + delay])
        clust_array = np.asarray(shape)
        dist = np.linalg.norm(clust_array - datum_array)
            
        #If the distance is the smallest so far, update.
        if dist < min_dist and np.max(np.asarray(datum)) == np.max(datum_array):
            min_dist = dist
    #Return the smallest distance
    return min_dist

if __name__ == "__main__":
    main()