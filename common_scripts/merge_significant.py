from scipy import shape
import numpy as np
import sys
import math
import os
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
If one shape is a shift of another, merge them.
"""
def main():

    file_path = sys.argv[1]
    output_path = sys.argv[2]
    merge_log = sys.argv[3]

    #Open the file containing the significant shapes.
    shape_list = []
    shape_anno = []
    shape_names = []
    file = open(file_path, 'r')
    mlog = open(merge_log, 'w')
    next_line = file.readline()
    while next_line:
        split_tabs = next_line.split("\t")
        shape_list.append([float(i) for i in split_tabs[2].split(",")])
        shape_anno.append(split_tabs[0])
        shape_names.append(split_tabs[1])
        next_line = file.readline()
            
    #Compare each shape to all shapes after it.
    #If the shapes should be merged and have the same annotation, shift the prior one as needed,
    #then delete the latter one.
    #If the shapes should be merged and one is unknown, keep the one that is not unknown.
    #If the shapes should be merged and have different annotations, remove both.
    #A threshold of 0.75 is used like in CoSBI
    length = len(shape_list)
    prev_j = 0
    for i in range(length - 1):
        j = i + 1
        while j < length:
            #If i and j meet the threshold, follow merging procedure.
            if should_merge(i, j, shape_list, 0.75):
                mlog.write("Match between " + shape_anno[i] + " and " + shape_anno[j])
                mlog.write("(" + shape_names[i] + ")" + " " + "(" + shape_names[j] + ")")
                mlog.write("\n")
                #If they have the same annotation, merge them, get the new length,
                #and leave the comparison index as is.
                if shape_anno[i] == shape_anno[j]:
                    shift_shape(i, j, shape_list, shape_names, shape_anno)
                    length = len(shape_list)
                #If one in unknown, keep the one that is not.
                #If they have different annotations, remove both to eliminate ambiguity.
                elif shape_anno[i] == "Unknown":
                    del shape_list[i]
                    del shape_anno[i]
                    del shape_names[i]
                    length = len(shape_list)
                elif shape_anno[j] == "Unknown":
                    del shape_list[j]
                    del shape_anno[j]
                    del shape_names[j]
                    length = len(shape_list)
                else:
                    del shape_list[j]
                    del shape_anno[j]
                    del shape_names[j]
                    del shape_list[i]
                    del shape_anno[i]
                    del shape_names[i]
                    length = len(shape_list)
            else:
                j += 1
                        
        #Print shifted shape shapes.
        outfile = open(output_path, 'w')
        for shape in range(len(shape_list)):
            outfile.write(shape_names[shape] + "\t" + shape_anno[shape] + "\t")
            for val in range(len(shape_list[shape])):
                outfile.write(str(shape_list[shape][val]))
                if val < len(shape_list[shape]) - 1:
                    outfile.write(",")
            outfile.write("\n")
        
    #Print message to user.
    print("Merging complete for all shapes in " + file_path)
    mlog.close()
    outfile.close()
    file.close()
    
#Determine whether two shapes should be merged.
def should_merge(shape1, shape2, shape_list, threshold):
    merge = False

    #Allow a shift of 1/4 of the shape dimension to each side.
    shift = 0.25 * len(shape_list[shape1])
    #Get the maximum cross-correlation metric.
    max_cc = 0
    for delay in range(int(-shift), int(shift)):
        try:
            cross_correlation = wsu.get_crosscorr(shape_list[shape1], shape_list[shape2], delay, threshold, 2, True, True, 0)
            if cross_correlation > max_cc:
                max_cc = cross_correlation

        #If the regions aren't valid for comparison at this delay, skip them.
        except ValueError:
            pass

    #If the maximum is above the threshold, then merge.
    if max_cc >= threshold:
        merge = True

    #Return decision about merging.
    return merge
    
#Use whichever shape has the most mass distributed at the center.
def shift_shape(shape, replacement, shape_list, names, annotations):
    
    #Determine which shape has its maximum closer to the center.
    max_index_shape = shape_list[shape].index(max(shape_list[shape]))
    max_index_replacement = shape_list[replacement].index(max(shape_list[replacement])) 
    
    #If the latter shape has its max at the center, replace the prior shape with it.
    if abs((len(shape_list[replacement]) / 2) - max_index_replacement) < abs((len(shape_list[shape]) / 2) - max_index_shape):
        shape_list[shape] = shape_list[replacement]
        
    #Remove the latter shape.
    del shape_list[replacement]
    del annotations[replacement]
    del names[replacement]

if __name__ == "__main__":
    main()