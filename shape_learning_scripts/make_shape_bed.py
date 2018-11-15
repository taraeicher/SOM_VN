#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

#Import numpy for calculations
import numpy as np

#Import sys for obtaining command line args.
import sys
sys.path.append(os.path.abspath("../annotation_scripts"))
import make_annotated_bed as mab

"""
Match each input to its closest shape. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    input = sys.argv[1]
    shape = sys.argv[2]
    output = sys.argv[3]
    minimum = float(sys.argv[4])

    #Create grids for labeling SOM nodes with shape indices.
    learned_shapes = []
    learned_shape_counts = []
    shape_names = []
    shape_count = 0
    
    #Open all input files.
    in_file = open(input, 'r')
    shape_file = open(shape, 'r')

    #Notify the user that files were found
    print("Using the input file:")
    print(in_file.name)

    #Notify the user that files were found
    print("Using the shape files:")
    print(shape_file.name)

    #Match the inputs to shapes and close the files.
    match_shapes(in_file.name, shape_file.name, output, minimum)
    in_file.close()
    shape_file.close()
    
"""
Match each input with the given window size to the nearest shape for that window size.
Print out the region with its corresponding shape to a BED file.
"""
def match_shapes(in_file_name, shape_file_name, out_dir, min_val):
    
    #Open input and shape files.
    in_file = open(in_file_name, "r")
    shape_file = open(shape_file_name, "r")
    
    #Open output file.
    out_file = open(out_dir, "w")
    
    #Read in shape data.
    shapes = []
    next_shape = shape_file.readline()
    counter = 0
    while next_shape:
        shapes.append([float(i) for i in next_shape.split(",")])
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
            
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [match, ambig, out_str] = match_region(input, shapes, min_val)
            
            #Print match to BED file. Format is:
            #chrom  start   end shape_num 1 - ambiguity
            #Score is the opposite of the ambiguity metric.
            out_file.write("chr" + labels[0] + "\t" + labels[1] + "\t" + labels[2] + "\t" + str(match) + "\t" + str(1 - ambig) + "\t" + out_str + "\n")

            #Read the next line in the file.            
            next_line = in_file.readline()
            
        #Print a message to the user.
        print("Files done")
    else:
        print("Only one shape. No annotation performed.")
    out_file.close()
    in_file.close()

"""
Find the closest match for the region in the list of shapes.
Return an ambiguity metric which measures how close the region
is to its nearest shape as opposed to other shapes.
"""
def match_region(region, shapes, min_val):

    #Create array to hold distances between region and shapes.
    max_crosscorr = 0
    match = 0
    opt_delay = 0
    crosscorr_list = []
    
    #For each shape, determine the region's distance from it.
    for i in range(len(shapes)):
        shape = shapes[i]
        crosscorr, d = mab.get_max_crosscorr(region, shape, min_val)
        crosscorr_list.append(crosscorr)
        if crosscorr_list[i] > max_crosscorr:
            max_crosscorr = crosscorr_list[i]
            match = i
            opt_delay = d

    #Calculate the ambiguity metric. Only do so if not all cross-correlations are zero.
    #If they are all zero, set ambiguity to 1 (completely ambiguous) and assign it to the
    #last shape (artificially created, with all 0.1).
    ambig = 1.0
    if not max_crosscorr == 0:
        ambig = mab.get_ambiguity(crosscorr_list, region)
    else:
        match = len(shapes) - 1
    
    #Write the shape values to a separate file.
    #This file will be used for shape validity analysis.
    comma = ","
    out_str = comma.join(str(e) for e in region[int(opt_delay):len(shapes[0]) + int(opt_delay)])
            
    #The ambiguity metric is the ratio of the closest shape's distance
    #to the distance of the second closest shape.
    return [match, ambig, out_str]
    
if __name__ == "__main__":
    main()