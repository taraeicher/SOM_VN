#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

#Import numpy for calculations
import numpy as np

#Import sys for obtaining command line args.
import sys
import os
sys.path.append(os.path.abspath("../annotation_scripts"))
import make_annotated_bed as mab

"""
Print the region and its magnitude to the nearest integer.
"""
def main():

    input = sys.argv[1]
    output = sys.argv[2]

    #Create grids for labeling SOM nodes with shape indices.
    learned_shapes = []
    learned_shape_counts = []
    shape_names = []
    shape_count = 0
    
    #Open all input files.
    in_file = open(input, 'r')

    #Notify the user that files were found
    print("Using the input file:")
    print(in_file.name)

    #Save the magnitude and close the files.
    match_shapes(in_file.name, output)
    in_file.close()
    
"""
Print out the region with its corresponding magnitude to a BED file.
"""
def match_shapes(in_file_name, out_dir):
    
    #Open input and shape files.
    in_file = open(in_file_name, "r")
    
    #Open output file.
    out_file = open(out_dir, "w")
    
    #Read in each line in the file and map it.
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
        magnitude = np.max(input)
        
        #Print match to BED file. Format is:
        #chrom  start   end shape_num 1 - ambiguity
        #Score is the opposite of the ambiguity metric.
        out_file.write("chr" + str(labels[0]) + "\t" + str(labels[1]) + "\t" + str(labels[2]) + "\t" + str(int(round(magnitude))) + "\n")

        #Read the next line in the file.            
        next_line = in_file.readline()
        
    #Print a message to the user.
    print("Files done")
    out_file.close()
    in_file.close()
    
if __name__ == "__main__":
    main()