#Import numpy for calculations
import numpy as np

#Import glob for obtaining the files.
import glob, os
import pickle as pkl

#Import sys for obtaining command line args.
import sys
from tqdm import tqdm
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
Match each input to its closest shape. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    input = sys.argv[1]
    shape = sys.argv[2]
    output_file = sys.argv[3]
    p_promoter = float(sys.argv[4])
    p_enhancer = float(sys.argv[5])
    p_repressor = float(sys.argv[6])

    #Create grids for labeling SOM nodes with shape indices.
    learned_shapes = []
    learned_shapes_counts = []
    learned_shapes.append(list())
    learned_shapes_counts.append(list())
    shape_names = [[]]
    shape_count = 0
    
    #Open all input files.
    in_file = pkl.load(open(input, 'rb'))
    shape_file = pkl.load(open(shape, 'rb'))

    #Match the inputs to shapes.
    match_shapes(in_file, shape_file, output_file, p_promoter, p_enhancer, p_repressor)

"""
Match each input to the nearest shape.
Print out the region with its corresponding shape to a BED file.
"""
def match_shapes(regions, shapes, out_file_name, p_promoter, p_enhancer, p_repressor):
        
    #Open output files.
    out_file = open(out_file_name, "w")
    
    #Read in each line in the file and map it.
    if(len(shapes) > 1):
        for region in tqdm(regions):
                   
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [match_label, score] = match_region(region.signals, shapes, p_promoter, p_enhancer, p_repressor)

            #Print match to BED file. Format is:
            #chrom  start   end shape_num 1 - ambiguity
            #Score is the opposite of the ambiguity metric.
            out_file.write("chr" + region.chromosome + "\t" + str(int(region.start)) + "\t" + str(int(region.end)) + "\t" + match_label + "\t" + str(score) + "\n")
            
        #Print a message to the user.
        print("Files done")
    else:
        print("Only one shape. No annotation performed.")
    out_file.close()

"""
Find the closest match for the input in the list of shapes.
Return an ambiguity metric which measures how close the input
is to its nearest shape as opposed to other shapes.
"""
def match_region(region, shapes, p_promoter, p_enhancer, p_repressor):

    #Create array to hold distances between region and shapes.
    max_crosscorr = 0
    match = 0
    opt_delay = 0
    crosscorr_list = []
    
    #For each shape, determine the region's distance from it.
    for i in range(len(shapes)):
        shape_assoc = shapes[i].shape.signals
        crosscorr, d = get_max_crosscorr(region, shape_assoc)
        crosscorr_list.append(crosscorr)
        if crosscorr_list[len(crosscorr_list) - 1] > max_crosscorr:
            max_crosscorr = crosscorr_list[len(crosscorr_list) - 1]
            match = i
            opt_delay = d
                
    # Get the annotation type comprising the maximum of regions
    # with this shape.
    label = None
    label = wsu.get_annotation(shapes[match], p_promoter, p_enhancer, p_repressor)
    percent_label = 0
    if label == "Promoter":
        percent_label = shapes[match].promoter_percentage
    elif label == "Enhancer":
        label = shapes[match].enhancer_percentage
    elif label == "Repressor":
        percent_label = shapes[match].repressed_percentage
    elif label == "Weak":
        percent_label = shapes[match].weak_percentage
        
    # Get the confidence of this annotation.
    return [label, percent_label * ((max_crosscorr + 1) / 2)]

"""
Find the minimum distance between the region and shifted shape.
"""
def get_max_crosscorr(region, shape):

    #Since the region is twice the size of the shape, search along the region and match at the best point.
    maximum = 0
    crosscorr = 0
    opt_delay = 0
    step = len(shape) / 10
    for i in range(0, 10):
        try:
            delay = i * step
            crosscorr = wsu.get_crosscorr(region, shape, int(delay))
  
            #If the distance is the smallest so far, update.
            if crosscorr > maximum:
                maximum = crosscorr
                opt_delay = delay
        except ValueError:
            pass
    #Return the smallest distance
    return [crosscorr, delay]

if __name__ == "__main__":
    main()