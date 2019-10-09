import sys
import plot_precision_recall as ppr
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
    our_bed = np.genfromtxt(sys.argv[1])
    pr_path = sys.argv[2]

    #Save precision and recall for the given experiment. 
    #Open the file to test that it exists.
    open_test = open(our_bed, "r")
    open_test.close()
    
    #Get all precision and recall values.
    [precision, recall] = get_all_precision_or_recall(our_bed)
    
    #Output reports for all types of annotations, including unknown.
    if len(precision) > 0:
        print_report(precision, recall, chrom, cell, win, pr_path)
    
def get_all_precision_or_recall(bed):

    annotations = ["Promoter", "Enhancer", "Repressed", "Weak"]
    length = len(annotations)

    #Get precision and recall for each type.
    [pred, gt] = get_labels_and_ground_truth(bed, annotations)
    precision = dict()
    recall = dict()
    if len(pred) > 0:
        for i in range(length):
            precision[i] = precision_score(gt[:, i], pred[:, i])
            recall[i] = recall_score(gt[:, i], pred[:, i])
        
    return [precision, recall]
    
#Get the percentage of the peak belonging to each ChromHMM annotation.
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
    
    #Do not move forward if the first line is blank in the sig file.
    try:
        sum_vec = np.zeros(4)
        
        #Keep track of regions with no ChromHMM annotations.
        #These regions will not be used in the analysis.
        not_annotated_count = 0
        count_in_region = 0
        for i in range(0, bed.shape[0]):            
                
            #Get the next element data.
            next_line = bed[i,:]
            current_start = int(next_line[1])
            current_end = int(next_line[2])
            a = next_line[8]
            anno_start = int(next_line[6])
            anno_end = int(next_line[7])
            our_anno = next_line[3]
            anno_length = int(next_line[9])
            
            #Get next signals if needed.
            #If we are still on the same region, don't get it.
            if current_start != prev_start:
                sum_vec = np.zeros(3)
                sig_i += 1
           
            #Add to the existing percentages.
            if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
                    sum_vec[0] += anno_length
            elif a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
                    sum_vec[1] += anno_length
            elif a == "13_ReprPC" or a == "14_ReprPCWk":
                    sum_vec[2] += anno_length
            elif a == "9_Het" or a == "15_Quies":
                    sum_vec[3] += anno_length
            
            #This case is when there is no annotation. Do not count it.
            else:
                not_annotated_count += 1
            count_in_region += 1
            
            #Add the ground truth and predicted value.
            next_start = current_start + 1
            if i + 1 < len(bed):
                next_start = int(bed[i + 1,:][1])
            if next_start != current_start and not_annotated_count != count_in_region:

                #Add another element to the ground truth and prediction vectors.
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
                #If it is unknown according to our analysis, do not consider it.
                #This includes cases where there is no annotation from ChromHMM or where
                #There is signal above the threshold but no ChromHMM annotation in the signal.
                if vec_pred[len(vec_pred) - 1][0] == 0 and vec_pred[len(vec_pred) - 1][1] == 0 and vec_pred[len(vec_pred) - 1][2] == 0:
                    del vec_pred[len(vec_pred) - 1]
                    del vec_gt[len(vec_pred) - 1]
                elif sum_vec[0] == 0 and sum_vec[1] == 0 and sum_vec[2] == 0:
                    del vec_pred[len(vec_pred) - 1]
                    del vec_gt[len(vec_pred) - 1]
                    
                #Set count and unannotated count to 0. Do the same for summation vec.
                not_annotated_count = 0
                count_in_region = 0
                
            #Get the previous data, if applicable.
            prev_start = current_start
            prev_end = current_end
            
        #Stack all values.
        final_stack_pred = np.stack(vec_pred)
        final_stack_gt = np.stack(vec_gt)
    except:
        pass
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
    
if __name__ == "__main__":
    main()