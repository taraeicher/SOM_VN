#Import numpy for calculations
import numpy as np

#Import glob for obtaining the files.
import glob, os

#Import sys for obtaining command line args.
import sys
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu
sys.path.append(os.path.abspath("../annotation_scripts"))
import make_annotated_bed as mab

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

    #Match the inputs to shapes and close the files.
    match_shapes_cutoff(in_file.name, shape_file.name, output_dir, wig_path, cutoff)
    in_file.close()
    shape_file.close()

"""
Match each input to the nearest shape.
Print out the region with its corresponding shape to a BED file.
"""
def match_shapes_cutoff(in_file_name, shape_file_name, out_dir, wig_name, cutoff):
    
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
            [match, ambig, crosscorr, out_str] = mab.match_region(input_scaled, shapes, shape_anno)
            
            #If the cross-correlation is within the threshold, assign the label as the correct label for the region.
            #Otherwise, assign it as unknown.
            anno_label = "Unknown"
            anno_name = "Unknown"
            if crosscorr >= cutoff and match != -1:
                anno_label = shape_anno[match]
                anno_name = shape_name[match]


            #Print match to BED file. Format is:
            #chrom  start   end cluster_num 1 - ambiguity
            #Score is the opposite of the ambiguity metric.
            if anno_label != "Unknown":
                out_file.write("chr" + labels[0] + "\t" + labels[1] + "\t" + labels[2] + "\t" + anno_label + "\t" + str(crosscorr) + "\n")

            #Read the next line in the file.            
            next_line = in_file.readline()
            
        #Print a message to the user.
        print("Files done")
    else:
        print("Only one cluster")
    out_file.close()
    in_file.close()

if __name__ == "__main__":
    main()