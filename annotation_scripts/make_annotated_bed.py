#Import numpy for calculations
import numpy as np

#Import glob for obtaining the files.
import glob, os
import pickle

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
    p_promoter = sys.argv[4]
    p_enhancer = sys.argv[5]
    p_repressor = sys.argv[6]

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

    #Notify the user that files were found
    print("Using the input file:")
    print(in_file.name)

    #Notify the user that files were found
    print("Using the shape file:")
    print(shape_file.name)

    #Match the inputs to shapes and close the files.
    match_shapes(in_file, shape_file, output_file, p_promoter, p_enhancer, p_weak, p_repressor)
    in_file.close()
    shape_file.close()

"""
Match each input to the nearest shape.
Print out the region with its corresponding shape to a BED file.
"""
def match_shapes(regions, out_file_name, shapes, p_promoter, p_enhancer, p_weak, p_repressor):
        
    #Open output files.
    out_file = open(out_file_name, "w")
    
    #Read in each line in the file and map it.
    if(len(shapes) > 1):
        for region in tqdm(regions):
                   
            #Match the data to the nearest shape and obtain the match and the ambiguity metric.
            [match, ambig, crosscorr, out_str] = match_region(region.signal, shapes, )
            
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

            #Read the next line in the file.            
            next_line = in_file.readline()
            
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
def match_region(region, shapes):

    #Create array to hold distances between region and shapes.
    max_crosscorr = 0
    match = 0
    opt_delay = 0
    crosscorr_list = []
    
    #For each shape, determine the region's distance from it.
    for i in range(len(shapes)):
        shape = shapes[i]
            crosscorr, d = get_max_crosscorr(region, shape)
            crosscorr_list.append(crosscorr)
            if crosscorr_list[len(crosscorr_list) - 1] > max_crosscorr:
                max_crosscorr = crosscorr_list[len(crosscorr_list) - 1]
                match = i
                opt_delay = d
                
    # Get the annotation type comprising the maximum of regions
    # with this shape.
    label = None
    argmax = wsu.get_annotation(match, p_promoter, p_enhancer, p_weak, p_repressor)
    if argmax == 0:
        label = "Promoter"
    elif argmax == 1:
        label = "Enhancer"
    elif argmax == 2:
        label = "Repressed"
    elif argmax == 3:
        label = "Weak"
        
    # Get the confidence of this annotation.
    total_anno = np.sum(match)
    percent_label = match[argmax] / total_anno
    return [label, percent_label]

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