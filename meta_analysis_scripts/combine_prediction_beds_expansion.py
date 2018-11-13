import numpy as np
import sys

BIN_SIZE = 50
"""
Combine the TSS file and our file to get a new file.
Follow "If this is a position-based promoter prediction OR
a shape-based promoter prediction, predict that it is a promoter.
Otherwise, leave it as is.
"""
def main():
    
    #Get list of clusters.
    tss_anno = np.genfromtxt(sys.argv[2], delimiter = "\t", dtype = str)
    our_anno = np.genfromtxt(sys.argv[1], delimiter = "\t", dtype = str)
    out_file = open(sys.argv[3], 'w')
    
    #Combine the annotations.
    make_new_file(tss_anno, our_anno, out_file)
    
    #Close the output file.
    out_file.close()

"""
For each line in the file, handle the following cases:
2. If the position-based prediction is not a
    promoter, output our region.
3. If the position-based prediction is a promoter, expand and output.
    a. If it can be expanded to fit between existing annotations,
        expand it. Position-based promoters must be expanded to 4000 bp. 
    b. Otherwise, replace the existing annotation.
Each time we output something, increment in that file.
"""    
def make_new_file(tss, ours, out):
    
    #Positions in the file
    chrom = 0
    start = 1
    end = 2
    annotation = 3
    
    #Counters for prediction-based and shape-based files.
    i = 0
    j = 0
    
    #Initialize start and end for TSS and our file.
    tss_line = tss[0]
    our_line = ours[0]
    previous_region = ""
    previous_promoter = ""
    
    #Loop through each file and consolidate.
    while i < len(tss) or j < len(ours):

        #If we've reached the end of our annotations, print and increment TSS.
        if j >= len(ours) and i < len(tss):
            if tss_line[annotation] == "Promoter":
                [new_line, overlaps_second] = get_expanded(tss_line, previous_region, "", previous_promoter, start, end, annotation)
                out.write(new_line + "\n")
                previous_promoter = new_line.split("\t")
            if i < len(tss):
                i += 1
        
        #If we've reached the end of the TSS annotations, print and increment the shape-based.
        elif i >= len(tss) and j < len(ours):
            out.write("\t".join(our_line) + "\n")
            previous_region = our_line
            if j < len(ours):
                j += 1
                
        #If the next shape-based annotation still comes before this TSS, increment shape-based only.
        elif int(our_line[end]) <= int(tss_line[start]) and tss_line[annotation] == "Promoter":
            out.write("\t".join(our_line) + "\n")
            previous_region = our_line
            j += 1
                
         #If the next shape-based annotation still comes before this TSS but TSS is not a promoter, increment both.
        elif int(our_line[end]) <= int(tss_line[start]) and tss_line[annotation] != "Promoter":
            if i < len(tss) - 1 and j < len(ours) - 1 and int(ours[j+1][end]) <= int(tss[i+1][start]):
                out.write("\t".join(our_line) + "\n")
                previous_region = our_line
                j += 1
            i += 1
            
        #If there is overlap or shape-based annotation comes last and TSS is not a promoter, increase TSS.
        elif int(our_line[end]) > int(tss_line[start]) and tss_line[annotation] != "Promoter":
            i += 1
                
        #If the next TSS comes before the next shape-based annotation, increment TSS only.
        #elif int(our_line[start]) > int(tss_line[end]) and tss_line[annotation] == "Promoter":
        elif tss_line[annotation] == "Promoter":
          
            #Skip everything that overlaps with the existing region.
            region_after = our_line
            while j < len(ours) and has_overlap(int(tss_line[start]), int(tss_line[end]), int(ours[j][start]), int(ours[j][end])):
                j += 1
                region_after = ""
                if j < len(ours):
                    our_line = ours[j]
                    region_after = our_line
            [new_line, overlaps_second] = get_expanded(tss_line, previous_region, region_after, previous_promoter, start, end, annotation)
            out.write(new_line + "\n")
            previous_promoter = new_line.split("\t")
            if overlaps_second:
                j += 1
            i += 1
                
        #If there is overlap or shape-based annotation comes last and TSS is a promoter, increment shape-based prediction.
        elif int(our_line[end]) > int(tss_line[start]) and tss_line[annotation] == "Promoter":
            j += 1
                
        #Get the next line.
        tss_line = tss[min(i, len(tss) - 1)]
        our_line = ours[min(j, len(ours) - 1)]
                   
"""
Determine whether there is overlap between two regions based on their start and end.
"""    
def has_overlap(start1, end1, start2, end2):              
    return_value = False
    if start1 >= start2 and start1 < end2:
        return_value = True
    elif start2 >= start1 and start2 < end1:
        return_value = True
    return return_value
   
"""
Find the next region that does not overlap with the current transcription start site.
"""    
# def find_next_region(tss, shape, tss_i, shape_i):

    # #As long as the regions still overlap, find the next shape-based annotation.
    # while has_overlap(int(tss[tss_i][start]), int(tss[tss_i][end]), int(shape[shape_i][start]), int(shape[shape_i][end])):
        # shape_i += 1
    # #Return the first shape-based annotation not to overlap.
    # return shape[shape_i]
    
"""
Expand the TSS promoter region so that it does not overlap existing regions.
If this cannot be done, send a notification to overwrite.
"""    
def get_expanded(line, before, after, prev_promoter, start, end, annotation):

    #Pad the promoter region with 500 bp on each side.
    padding = 500
    expanded_size = 4000
    new_start = max(0, int(line[start]) - padding)
    if "\t".join(prev_promoter) != "":
        new_start = max(new_start, int(prev_promoter[start]))
    new_end = new_start + expanded_size
    new_line = line
     
    #Shift backward or forward until an optimal region is found.
    does_overlap_prev = False
    does_overlap_promoter = False
    if "\t".join(before) != "":
        does_overlap_prev = has_overlap(new_start, new_end, int(before[start]), int(before[end]))
    if "\t".join(prev_promoter) != "":
        does_overlap_promoter = has_overlap(new_start, new_end, int(prev_promoter[start]), int(prev_promoter[end])) and prev_promoter[annotation] == "Promoter"
        
    #If we overlap a previous promoter, append to the end of that promoter.
    if does_overlap_promoter:
        new_start = max(0, int(line[start]) - padding, int(prev_promoter[end]))
        new_end = new_start + expanded_size
        
    #If we overlap the first region after expanding, shift until we don't.
    elif does_overlap_prev:
        while does_overlap_prev and new_end - BIN_SIZE >= int(line[end]):
            new_start += BIN_SIZE
            new_end += BIN_SIZE
            does_overlap_prev = has_overlap(new_start, new_end, int(before[start]), int(before[end]))
            
    #Get overlap info for region after.
    does_overlap_post = False
    if "\t".join(after) != "":
        does_overlap_post = has_overlap(new_start, new_end, int(after[start]), int(after[end]))
                
    #Return the optimal region.
    return ["\t".join([line[0], str(new_start), str(new_end), line[3], str(1.0)]), does_overlap_post]
    
if __name__ == "__main__":
    main()