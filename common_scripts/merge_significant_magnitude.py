import numpy as np
import sys
import math
import os
import glob
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
If one mag is a shift of another, merge them.
"""
def main():

    file_path = sys.argv[1]
    output_path = sys.argv[2]
    merge_log = sys.argv[3]

    #Open the file containing the significant mags.
    mag_list = []
    mag_anno = []
    mag_names = []
    mlog = open(merge_log, 'w')
    
    for f in glob.glob(file_path + "/*"):
        file = open(f, 'r')
        next_line = file.readline()
        while next_line:
            split_tabs = next_line.split("\t")
            mag_list.append(int(split_tabs[2]))
            mag_anno.append(split_tabs[0])
            mag_names.append(split_tabs[1])
            next_line = file.readline()
        file.close()
       
    #Compare each mag to all mags after it.
    #If the mags should be merged and have the same annotation, shift the prior one as needed,
    #then delete the latter one.
    #If the mags should be merged and one is unknown, keep the one that is not unknown.
    #If the mags should be merged and have different annotations, remove both.
    #A threshold of 0.75 is used like in CoSBI
    outfile = None
    length = len(mag_list)
    prev_j = 0
    i = 0
    while i in range(length - 1):
        j = i + 1
        increase_i = True
        while j < length:
            #If i and j meet the threshold, follow merging procedure.
            if mag_list[i] == mag_list[j]:
                mlog.write("Match between " + mag_anno[i] + " and " + mag_anno[j])
                mlog.write("(" + mag_names[i] + ")" + " " + "(" + mag_names[j] + ")")
                mlog.write("\n")
                #If they have the same annotation, merge them, get the new length,
                #and leave the comparison index as is.
                if mag_anno[i] == mag_anno[j]:
                    del mag_list[j]
                    length = len(mag_list)
                #If one in unknown, keep the one that is not.
                #If they have different annotations, remove both to eliminate ambiguity.
                elif mag_anno[i] == "Unknown":
                    del mag_list[i]
                    del mag_anno[i]
                    del mag_names[i]
                    increase_i = False
                    length = len(mag_list)
                elif mag_anno[j] == "Unknown":
                    del mag_list[j]
                    del mag_anno[j]
                    del mag_names[j]
                    length = len(mag_list)
                else:
                    del mag_list[j]
                    del mag_anno[j]
                    del mag_names[j]
                    del mag_list[i]
                    del mag_anno[i]
                    del mag_names[i]
                    increase_i = False
                    length = len(mag_list)
            else:
                j += 1
                
        if increase_i == True:
            i += 1
        
    #Print shifted mag mags.
    outfile = open(output_path, 'w')
    for mag in range(len(mag_list)):
        outfile.write(mag_names[mag] + "\t" + mag_anno[mag] + "\t")
        outfile.write(str(mag_list[mag]))
        outfile.write("\n")
        
    #Print message to user.
    print("Merging complete for all mags in " + file_path)
    mlog.close()
    outfile.close()

if __name__ == "__main__":
    main()