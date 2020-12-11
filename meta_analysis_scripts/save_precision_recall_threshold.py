import sys
import plot_precision_recall_nobaselines as ppr
import numpy as np
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
    pr_path = sys.argv[4]
    all_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
    precision_recall_list = []

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
            [chrom, precision, recall] = get_all_precision_or_recall(our_bed, our_sig, wig, chrom)
            precision_recall_list.append([chrom, precision, recall])
            
        except OSError:
            pass
    print_report(precision_recall_list, pr_path)
    
#Plot the ROC curve based on ground truth and prediction for each cutoff. Plot separate lines for each
#annotation type. Consolidate all chromosomes.
#def get_all_precision_or_recall(bed, sig, tss_bed, tss_sig, or_bed, or_sig, perm_bed, perm_sig, rpkm_bed, rpkm_sig, wig, chrom, win, cell):
def get_all_precision_or_recall(bed, sig, wig, chrom):

    #Get actual annotation and ground truth for all annotations and for all unannotated regions.
    threshold = 5
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
        
    return [chrom, precision, recall]
    
"""
Print precision and recall for all chromosomes.
""" 
def print_report(precision_recall_list, pr):

    report = open(pr, "w")
    
    #All shape-based predictions
    report.write("Chrom,Promoter_Precision,Promoter_Recall,Enhancer_Precision,Enhancer_Recall,Weak_Precision,Weak_Recall\n")
    for j in range(len(precision_recall_list)):
        report.write(precision_recall_list[j][0] + ",")
        report.write(str(precision_recall_list[j][1][0]) + "," + str(precision_recall_list[j][2][0]) + ",")
        report.write(str(precision_recall_list[j][1][1]) + "," + str(precision_recall_list[j][2][1]) + ",")
        report.write(str(precision_recall_list[j][1][2]) + "," + str(precision_recall_list[j][2][2]) + "\n")
    
    #Close the report.
    report.close()
    
if __name__ == "__main__":
    main()