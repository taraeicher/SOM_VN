import numpy as np
import scipy as sp
import sys
import math
import common_ops as ops
from scipy import stats
from scipy import interp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import common_ops as ops
import glob, os
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve


"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the chromHMM and annotated files.
    bed = sys.argv[1]
    tss_bed = sys.argv[2]
    roc_path = sys.argv[3]
    pr_path = sys.argv[4]
    unknown_path = sys.argv[5]
    wig = sys.argv[6]
    sig = sys.argv[7]
    tss_sig = sys.argv[8]
    chrom = sys.argv[9]
    win = sys.argv[10]
    cell = sys.argv[11]

    #Plot ROC curve for each chromosome. Plot curve for each annotation separately.
    output_reports(bed, sig, wig, tss_bed, tss_sig, roc_path, unknown_path, pr_path, chrom, win, cell)
   
#Plot the ROC curve based on ground truth and prediction for each cutoff. Plot separate lines for each
#annotation type. Consolidate all chromosomes.
def output_reports(bed, sig, wig, tss_bed, tss_sig, roc_path, unknown_path, report_path, chrom, win, cell):

    #Get actual annotation and ground truth for all annotations and for all unannotated regions.
    threshold = ops.get_intensity_percentile(0.75, open(wig, 'r'))
    annotations = ["Promoter", "Enhancer", "Weak", "Polycomb", "Other"]
    length = len(annotations)
    ground_truth_list = []
    predicted_list = []
    [pred, gt] = get_labels_and_ground_truth(bed, sig, wig, annotations, threshold)  

    #Get ROC curve for each annotation type.
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    precision = dict()
    recall = dict()
    for i in range(length):
        fpr[i], tpr[i], _ = roc_curve(gt[:, i], pred[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
        precision[i], recall[i], _ = precision_recall_curve(gt[:, i], pred[:, i])
        
    #Get ROC curve for TSS-based promoters.
    [pred_promoter, gt_promoter] = get_tss_labels_and_ground_truth(tss_bed, tss_sig, wig, ["Promoter", "Not_Promoter"], threshold) 
    fpr["tss"], tpr["tss"], _ = roc_curve(gt_promoter[:,0], pred_promoter[:,0])
    roc_auc["tss"] = auc(fpr["tss"], tpr["tss"])
    precision["tss"], recall["tss"], _ = precision_recall_curve(gt_promoter[:, 0], pred_promoter[:, 0])
    
    #Plot all ROC curves.
    plt.figure()
    colors = ['darkgreen', 'orangered', 'darkred', 'blue', 'slategray']
    for i, color in zip(range(length), colors):
        plt.plot(fpr[i], tpr[i], color=color, linewidth=4,
             label=annotations[i] + ' ROC curve(area = {0:0.2f})'
               ''.format(roc_auc[i]))
    plt.plot(fpr["tss"], tpr["tss"], color = 'black', linewidth = 2, label = "TSS-Based Promoter")
    plt.plot([0, 1], [0, 1], 'k--', linewidth=4)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Curves for Window Size ' + win + ' - ' + cell + " Chromosome " + chrom)
    plt.legend(loc="lower right")
    plt.savefig(roc_path + ".png")
    
    #Print report. Give the details of the data set and the TP/FP for:
    #Promoter
    #Enhancer
    #Weak
    #Polycomb
    #Other
    #Promoter TSS
    report = open(report_path, "w")
    report.write(cell + "\n")
    report.write(chrom + "\n")
    report.write(win + "\n")
    for i in range(length):
        report.write(str(precision[i][1]) + "\t" + str(recall[i][1]) + "\n")
    report.write(str(precision["tss"][1]) + "\t" + str(recall["tss"][1]) + "\n")
    report.close()
    
    #Print plots of unknown percentages for each chromosome.
    print_unknown_percentages(bed, sig, wig, unknown_path, pred.shape[0], chrom, win, cell, threshold)

#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def get_labels_and_ground_truth(bed_file, sig_file, wig, annotations, threshold):
        
    #Set up percentage matrix.
    vec_pred = list()
    vec_gt = list()
    
    #Get scores and labels for each bed file.
    bed = np.genfromtxt(bed_file , delimiter='\t', dtype = str)
    sigs = np.genfromtxt(sig_file, delimiter = ',', dtype = float)

    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    prev_start = -1
    prev_end = -1
    sig = sigs[0,:]
    sig_i = -1
    sum_vec = np.zeros(5)
    
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
            sig_i += 1
            sig = sigs[sig_i,:]
            sum_vec = np.zeros(5)
       
        #Add to the existing percentages.
        #If the region has peaks, consider only regions above peak threshold.
        #If no peaks exist, consider entire region.
        total_peak_size = count_above(threshold, "", sig, current_start, current_end, current_start, current_end)
        if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
            if total_peak_size > 0:
                sum_vec[0] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[0] += anno_length
        elif a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
            if total_peak_size > 0:
                sum_vec[1] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[1] += anno_length
        elif a == "13_ReprPC" or a == "ReprPCWk":
            if total_peak_size > 0:
                sum_vec[3] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[3] += anno_length
        elif a == "9_Het" or a == "15_Quies":
            if total_peak_size > 0:
                sum_vec[2] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[2] += anno_length
        elif anno_length != 0:
            if total_peak_size > 0:
                sum_vec[4] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[4] += anno_length
        #This case is when there is no annotation. Do not count it.
        else:
            not_annotated_count += 1
        count_in_region += 1
        
        #Add the ground truth and predicted value based on the region with the maximum count above the threshold.
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
            if vec_pred[len(vec_pred) - 1][0] == 0 and vec_pred[len(vec_pred) - 1][1] == 0 and vec_pred[len(vec_pred) - 1][2] == 0 and vec_pred[len(vec_pred) - 1][3] == 0 and vec_pred[len(vec_pred) - 1][4] == 0:
                del vec_pred[len(vec_pred) - 1]
                del vec_gt[len(vec_pred) - 1]
            elif sum_vec[0] == 0 and sum_vec[1] == 0 and sum_vec[2] == 0 and sum_vec[3] == 0 and sum_vec[4] == 0:
                del vec_pred[len(vec_pred) - 1]
                del vec_gt[len(vec_pred) - 1]
                
            #Set count and unannotated count to 0. Do the same for summation vec.
            not_annotated_count = 0
            count_in_region = 0
            
        #Get the previous data, if applicable.
        prev_start = current_start
        prev_end = current_end
    #Return value.
    return [np.stack(vec_pred), np.stack(vec_gt)]
    
#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def get_tss_labels_and_ground_truth(bed_file, sig_file, wig, annotations, threshold):
    
    #Set up counts of promoter and non-promoter.
    vec_pred = list()
    vec_gt = list()
    
    bed = np.genfromtxt(bed_file , delimiter='\t', dtype = str)
    sigs = np.genfromtxt(sig_file, delimiter = ',', dtype = float)
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    prev_start = -1
    prev_end = -1
    sig = sigs[0,:]
    sig_i = -1
    sum_vec = np.zeros(2)
    
    #Keep track of regions with no ChromHMM annotations.
    #These regions will not be used in the analysis.
    not_annotated_count = 0
    count_in_region = 0
    
    for i in range(0, bed.shape[0]):            
            
        #Get the next element data.
        next_line = bed[i,:]
        current_start = int(next_line[1])
        current_end = int(next_line[2])
        a = next_line[7]
        anno_start = int(next_line[5])
        anno_end = int(next_line[6])
        our_anno = next_line[3]
        anno_length = int(next_line[8])
        
        #Get next signals if needed.
        #If we are still on the same region, don't get it.
        if current_start != prev_start:
            sig_i += 1
            sig = sigs[sig_i,:]
            sum_vec = np.zeros(2)
            
        #Add to the existing percentages.
        #If the region has peaks, consider only regions above peak threshold.
        #If no peaks exist, consider entire region.
        total_peak_size = count_above(threshold, "", sig, current_start, current_end, current_start, current_end)
        if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
            if total_peak_size > 0:
                sum_vec[0] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[0] += anno_length
        elif anno_length != 0:
            if total_peak_size > 0:
                sum_vec[1] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[1] += anno_length
        #This case is when there is no annotation. Do not count it.
        else:
            not_annotated_count += 1
        count_in_region += 1
          
        #Add the ground truth and predicted value based on the region with the maximum count above the threshold.
        next_start = current_start + 1
        if i + 1 < len(bed):
            next_start = int(bed[i + 1,:][1])
        if next_start != current_start and not_annotated_count != count_in_region:
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
            if sum_vec[0] == 0 and sum_vec[1] == 0:
                del vec_pred[len(vec_pred) - 1]
                del vec_gt[len(vec_pred) - 1]
            #Set count and unannotated count to 0.
            not_annotated_count = 0
            count_in_region = 0
        #Get the previous data, if applicable.
        prev_start = current_start
        prev_end = current_start
    #Return value.
    return [np.stack(vec_pred), np.stack(vec_gt)]
    
#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def print_unknown_percentages(bed_file, sig_file, wig, out_path, total_regions, chrom, win, cell, threshold):
    
    #Set up counts of promoter and non-promoter.
    vec_gt = list()
  
    bed = np.genfromtxt(bed_file , delimiter='\t', dtype = str)
    sigs = np.genfromtxt(sig_file, delimiter = ',', dtype = float)
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    prev_start = -1
    prev_end = -1
    sig = sigs[0,:]
    sig_i = -1
    sum_vec = np.zeros(2)
    
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
            sig_i += 1
            sig = sigs[sig_i,:]
            sum_vec = np.zeros(5)
            
        #If this is an unknown region, see which ChromHMM annotation makes up the bulk of the region.
        if our_anno == "Unknown":
            #Add to the existing percentages.
            #If the region has peaks, consider only regions above peak threshold.
            #If no peaks exist, consider entire region.
            total_peak_size = count_above(threshold, "", sig, current_start, current_end, current_start, current_end)
            if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
                if total_peak_size > 0:
                    sum_vec[0] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[0] += anno_length
            elif a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
                if total_peak_size > 0:
                    sum_vec[1] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[1] += anno_length
            elif a == "13_ReprPC" or a == "ReprPCWk":
                if total_peak_size > 0:
                    sum_vec[3] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[3] += anno_length
            elif a == "9_Het" or a == "15_Quies":
                if total_peak_size > 0:
                    sum_vec[2] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[2] += anno_length
            elif anno_length != 0:
                if total_peak_size > 0:
                    sum_vec[4] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[4] += anno_length
            #This case is when there is no annotation. Do not count it.
            else:
                not_annotated_count += 1
            count_in_region += 1
              
            #Add the ground truth and predicted value based on the region with the maximum count above the threshold.
            next_start = current_start + 1
            if i + 1 < len(bed):
                next_start = int(bed[i + 1,:][1])
            if next_start != current_start and not_annotated_count != count_in_region:
                vec_gt.append(np.zeros(len(sum_vec)))
                max = np.argmax(sum_vec)
                for sum in range(0, len(sum_vec)):
                    #Add ground truth.
                    if sum == max:
                        vec_gt[len(vec_gt) - 1][sum] = 1
                    else:
                        vec_gt[len(vec_gt) - 1][sum] = 0
                if sum_vec[0] == 0 and sum_vec[1] == 0 and sum_vec[2] == 0 and sum_vec[3] == 0 and sum_vec[4] == 0:
                    del vec_gt[len(vec_gt) - 1]
                #Set count and unannotated count to 0.
                not_annotated_count = 0
                count_in_region = 0

        #Get the previous data, if applicable.
        prev_start = current_start
        prev_end = current_start
            
    #Write a report with the percentage of unknowns and, of the unknowns, the percentage of each annotation type.
    out = open(out_path, "w")
    ground_truth = np.stack(vec_gt)
    unknown_percent = ground_truth.shape[0] / sig_i
    count_all = np.sum(ground_truth, axis = 0)
    out.write(str(unknown_percent) + "\n")
    percent_promoters = count_all[0] / ground_truth.shape[0]
    out.write(str(percent_promoters) + "\n")
    percent_enhancers = count_all[1] / ground_truth.shape[0]
    out.write(str(percent_enhancers) + "\n")
    percent_weak = count_all[2] / ground_truth.shape[0]
    out.write(str(percent_weak) + "\n")
    percent_polycomb = count_all[3] / ground_truth.shape[0]
    out.write(str(percent_polycomb) + "\n")
    percent_other = count_all[4] / ground_truth.shape[0]
    out.write(str(percent_other) + "\n")
 
def count_above(threshold, annotation, signal, start, end, start_anno, end_anno):
    count = 0
    start_idx = start
    for sig in signal:
        is_between_anno = (start_anno <= start_idx) and (start_idx <= end_anno)
        if sig > threshold and (is_between_anno or annotation == ""):
            count += BIN_SIZE
        start_idx += BIN_SIZE
    return count
    
if __name__ == "__main__":
    main()