import sys
import common_ops as ops
import plot_precision_recall as ppr
import numpy as np
import glob, os
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the chromHMM and annotated files.
    pr_path = sys.argv[4]
    win = sys.argv[5]
    cell = sys.argv[6]
    all_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    #Plot ROC curve for each chromosome. Plot curve for each annotation separately.
    c = 0
    for chrom in all_chroms:
        #Get files for each category.
        our_bed = sys.argv[1] + "anno" + str(chrom) + ".bed"
        our_sig = sys.argv[2] + "clusters_anno" + str(chrom)
        wig = sys.argv[3] + str(chrom) + ".wig"
        try:     
            #Open the file to test that it exists.
            open_test = open(our_bed, "r")
            open_test.close()
            
            #Get all precision and recall values.
            [precision, recall] = get_all_precision_or_recall(our_bed, our_sig, wig, chrom, win, cell)
            
            #Output reports for all types of annotations, including unknown.
            if len(precision) > 0:
                print_report(precision, recall, chrom, cell, win, pr_path)

        except FileNotFoundError:
            pass
    
#Plot the ROC curve based on ground truth and prediction for each cutoff. Plot separate lines for each
#annotation type. Consolidate all chromosomes.
#def get_all_precision_or_recall(bed, sig, tss_bed, tss_sig, or_bed, or_sig, perm_bed, perm_sig, rpkm_bed, rpkm_sig, wig, chrom, win, cell):
def get_all_precision_or_recall(bed, sig, wig, chrom, win, cell):

    #Get actual annotation and ground truth for all annotations and for all unannotated regions.
    threshold = ops.get_intensity_percentile(0.75, open(wig, 'r'))
    annotations = ["Promoter", "Enhancer", "Weak"]
    length = len(annotations)

    #Get precision and recall for each type.
    [pred, gt] = ppr.get_labels_and_ground_truth(bed, sig, wig, annotations, threshold)
    precision = dict()
    recall = dict()
    if len(pred) > 0:
        for i in range(length):
            precision[i] = precision_score(gt[:, i], pred[:, i])
            recall[i] = recall_score(gt[:, i], pred[:, i])
        
    return [precision, recall]
    
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
    
if __name__ == "__main__":
    main()