"""
This script computes the precision and recall values for each
peak in a BED file. Precision and recall are computed separately
for each RE.

Required parameters:
1. An intersected BED file with the following format:
    peak_start  peak_end    peak_annotation chromhmm_start  chromhmm_end    chromhmm_annotation overlap_length
2. The file where you wish to store the precision and recall in CSV format.
"""

import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import glob, os
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the chromHMM and annotated files.
    our_bed = pd.read_csv(sys.argv[1], sep = "\t")
    pr_path = sys.argv[2]
    
    #Get all precision and recall values.
    get_all_precision_and_recall(our_bed).to_csv(open(pr_path, "w"))
    
"""
Get the precision and recall for each RE.
"""
def get_all_precision_and_recall(bed):

    annotations = ["Promoter", "Enhancer", "Repressor", "Weak"]
    length = len(annotations)

    #Get precision and recall for each type.
    [pred, gt] = get_labels_and_ground_truth(bed, annotations)
    promoter_pr = [precision_score(gt[:,0], pred[:,0]), recall_score(gt[:,0], pred[:,0])]
    enhancer_pr = [precision_score(gt[:,1], pred[:,1]), recall_score(gt[:,1], pred[:,1])]
    repressor_pr = [precision_score(gt[:,2], pred[:,2]), recall_score(gt[:,2], pred[:,2])]
    weak_pr = [precision_score(gt[:,3], pred[:,3]), recall_score(gt[:,3], pred[:,3])]
    pr_all = pd.DataFrame({"Promoter": promoter_pr,
        "Enhancer": enhancer_pr,
        "Repressor": repressor_pr,
        "Weak": weak_pr},
        index = ["Precision", "Recall"],
        columns = annotations)
        
    return pr_all
    
"""
Get the percentage of the chromosome belonging to each ChromHMM annotation.
"""
def get_labels_and_ground_truth(bed, annotations):
    
    #Set up percentage matrix.
    vec_pred = list()
    vec_gt = list()
    final_stack_pred = np.empty((0, 0))
    final_stack_gt = np.empty((0, 0))

    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    prev_start = -1
    prev_end = -1
    sum_vec = np.zeros(4)
    
    #Keep track of regions with no ChromHMM annotations.
    #These regions will not be used in the analysis.
    not_annotated_count = 0
    count_in_region = 0
    for i in tqdm(range(0, bed.shape[0])):            
            
        #Get the next peak-ChromHMM overlap.
        current_start = bed.iloc[i,1]
        current_end = bed.iloc[i,2]
        a = bed.iloc[i,8]
        anno_start = bed.iloc[i,6]
        anno_end = bed.iloc[i,7]
        our_anno = bed.iloc[i,3]
        anno_length = bed.iloc[i,9]
        
        #If this ChromHMM annotation is the first to overlap
        #with this peak, zero out the sum vector and start summing
        #ChromHMM annotation size counts for this peak.
        if current_start != prev_start:
            sum_vec = np.zeros(4)
       
        #Add to the existing peak counts of ChromHMM annotations for that peak.
        [sum_vec, not_annotated_count, count_in_region] = add_chromhmm_annotation(sum_vec, anno_length, a, not_annotated_count, count_in_region)
        
        #If we are done with this peak, compute ground truth vs predicted.
        next_start = current_start + 1
        if i + 1 < len(bed):
            next_start = int(bed.iloc[i + 1,1])
        if next_start != current_start and not_annotated_count != count_in_region:

            #Add another element to the ground truth and prediction vectors.
            [vec_gt, vec_pred] = add_another_gt_pred(vec_gt, vec_pred, sum_vec, annotations, our_anno)
                
            #Set count and unannotated count to 0. Do the same for summation vec.
            not_annotated_count = 0
            count_in_region = 0
            
        #Get the previous data, if applicable.
        prev_start = current_start
        prev_end = current_end
            
        #Stack all values.
        final_stack_pred = np.stack(vec_pred)
        final_stack_gt = np.stack(vec_gt)
        
    #Return value.
    return [final_stack_pred, final_stack_gt]
    
"""
Print precision and recall for all chromosomes.
""" 
def print_report(precision, recall, chrom, win, cell, pr):
    #Print the cell line, chromosome, and window information
    report = open(pr + "/report_" + chrom, "w")
    report.write(cell + "\n")
    report.write(chrom + "\n")
    report.write(win + "\n\n")
    
    #All shape-based predictions
    for i in range(0, 3):
        report.write(str(precision[i]) + "\t" + str(recall[i]) + "\n")
    report.write("\n")
    
    #Close the report.
    report.close()
    
"""
Using the next known overlap between peak and ChromHMM, add to the summation vector of bp with
each mnemonic in this peak and track the total number of overlaps and overlaps without mnemonics
in this peak.
"""
def add_chromhmm_annotation(sum_vec, overlap_length, chromhmm_annotation, not_annotated_count, count_in_region):
    promoters = ["1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"]
    enhancers = ["6_EnhG", "7_Enh", "12_EnhBiv"]
    repressors = ["13_ReprPC", "14_ReprPCWk"]
    weak = ["9_Het", "15_Quies"]
    if chromhmm_annotation in promoters:
        sum_vec[0] += overlap_length
    elif chromhmm_annotation in enhancers:
        sum_vec[1] += overlap_length
    elif chromhmm_annotation in repressors:
        sum_vec[2] += overlap_length
    elif chromhmm_annotation in weak:
        sum_vec[3] += overlap_length
    #This case is when there is no annotation. Do not count it.
    else:
        not_annotated_count += 1
    count_in_region += 1
    return [sum_vec, not_annotated_count, count_in_region]
    
"""
Add the ground truth annotation of the peak as the
annotation that is most prominent inside the peak.
Add the annotation predicted using our method to the
vector of predictions.
"""
def add_another_gt_pred(vec_gt, vec_pred, sum_vec, annotations, our_anno):
    vec_gt.append(np.zeros(len(sum_vec)))
    vec_pred.append(np.zeros(len(sum_vec)))
    max = np.argmax(sum_vec)
    for sum in range(0, len(sum_vec)):
        #Add ground truth.
        if sum == max:
            vec_gt[len(vec_gt) - 1][sum] = 1
        else:
            vec_gt[len(vec_gt) - 1][sum] = 0 
        #Add predictions.
        if annotations[sum] == our_anno:
            vec_pred[len(vec_pred) - 1][sum] = 1
        else:
            vec_pred[len(vec_pred) - 1][sum] = 0 
    return [vec_gt, vec_pred]
    
if __name__ == "__main__":
    main()