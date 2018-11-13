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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import common_ops as ops
import glob, os
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_recall_curve
import seaborn as sns

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the chromHMM and annotated files.
    plot_out = sys.argv[11]
    pr_path = sys.argv[12]
    unknown_path = sys.argv[13]
    win = sys.argv[15]
    cell = sys.argv[16]
    src = sys.argv[17]
    #indices_to_highlight = [3, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    indices_to_highlight = []
    all_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    precision_all = np.zeros((3, len(all_chroms)))
    recall_all = np.zeros((3, len(all_chroms)))
    #precision_tss = np.zeros((2, len(all_chroms)))
    #recall_tss = np.zeros((2, len(all_chroms)))
    #precision_or_all = np.zeros((3, len(all_chroms)))
    #recall_or_all = np.zeros((3, len(all_chroms)))
    precision_and_all = np.zeros((3, len(all_chroms)))
    recall_and_all = np.zeros((3, len(all_chroms)))
    #precision_perm_all = np.zeros((3, len(all_chroms)))
    #recall_perm_all = np.zeros((3, len(all_chroms)))
    #precision_rpkm_all = np.zeros((3, len(all_chroms)))
    #recall_rpkm_all = np.zeros((3, len(all_chroms)))
    predictions_all = []
    ground_truth_all = []

    #Plot ROC curve for each chromosome. Plot curve for each annotation separately.
    c = 0
    for chrom in all_chroms:
        #Get files for each category.
        our_bed = sys.argv[1] + "anno" + str(chrom) + ".bed"
        our_sig = sys.argv[5] + "clusters_anno" + str(chrom)
        tss_bed = sys.argv[2] + "anno" + str(chrom) + ".bed"
        tss_sig = sys.argv[6] + "clusters" + str(chrom)
        or_bed = sys.argv[1] + "anno" + str(chrom) + "_tss.bed"
        or_sig = sys.argv[7] + "clusters" + str(chrom)
        and_bed = sys.argv[1] + "anno" + str(chrom) + "_tss_and.bed"
        and_sig = sys.argv[8] + "clusters_anno" + str(chrom)
        perm_bed = sys.argv[4] + "anno" + str(chrom) + ".bed"
        perm_sig = sys.argv[10] + "clusters_anno" + str(chrom)
        rpkm_bed = sys.argv[3] + "anno" + str(chrom) + "_final.bed"
        rpkm_sig = sys.argv[9] + "clusters_anno" + str(chrom)
        wig = sys.argv[14] + str(chrom) + ".wig"
        
        #Get all precision and recall values.
        #[precision, recall, precision_or, recall_or, precision_perm, recall_perm, precision_rpkm, recall_rpkm, total, threshold, predictions, ground_truth] = get_all_precision_or_recall(our_bed, our_sig, tss_bed, tss_sig, or_bed, or_sig, perm_bed, perm_sig, rpkm_bed, rpkm_sig, wig, chrom, win, cell)
        [precision, recall, precision_and, recall_and, total, threshold, predictions, ground_truth] = get_all_precision_or_recall(our_bed, our_sig, and_bed, and_sig, wig, chrom, win, cell)
        predictions_all.append(predictions)
        ground_truth_all.append(np.asarray(ground_truth))
        
        for i in range(0,3):
            precision_all[i,c] = precision[i]
            recall_all[i,c] = recall[i]
            #precision_or_all[i,c] = precision_or[i]
            #recall_or_all[i,c] = recall_or[i]
            precision_and_all[i,c] = precision_and[i]
            recall_and_all[i,c] = recall_and[i]
            #precision_perm_all[i,c] = precision_perm[i]
            #recall_perm_all[i,c] = recall_perm[i]
            #precision_rpkm_all[i,c] = precision_rpkm[i]
            #recall_rpkm_all[i,c] = recall_rpkm[i]
            
        # for i in range(0,2):
            # precision_tss[i,c] = precision["tss"]
            # recall_tss[i,c] = recall["tss"]
        c += 1
        
        #Output reports for all types of annotations, including unknown.
        #print_report(precision, recall, precision_or, recall_or, precision_perm, recall_perm, precision_rpkm, recall_rpkm, chrom, cell, win, pr_path)
        
        #Print plots of unknown percentages for each chromosome.
        #print_unknown_percentages(our_bed, our_sig, wig, unknown_path, total, chrom, win, cell, threshold)
        
    #Save a scatterplot with all precision and recall values.
    #save_scatterplot(precision_all, recall_all, precision_tss, recall_tss, precision_or_all, recall_or_all, precision_and_all, recall_and_all, precision_rpkm_all, recall_rpkm_all, plot_out, cell, src, indices_to_highlight, all_chroms, True)
    save_scatterplot(precision_all, recall_all, precision_and_all, recall_and_all, plot_out, cell, src, indices_to_highlight, all_chroms, True)
    
    #Save a heatmap of all mispredictions.
    save_misprediction_heatmap(np.vstack(predictions_all[:22]), np.vstack(ground_truth_all[:22]), plot_out, cell, src)
    
#Plot the ROC curve based on ground truth and prediction for each cutoff. Plot separate lines for each
#annotation type. Consolidate all chromosomes.
#def get_all_precision_or_recall(bed, sig, tss_bed, tss_sig, or_bed, or_sig, perm_bed, perm_sig, rpkm_bed, rpkm_sig, wig, chrom, win, cell):
def get_all_precision_or_recall(bed, sig, and_bed, and_sig, wig, chrom, win, cell):

    #Get actual annotation and ground truth for all annotations and for all unannotated regions.
    threshold = ops.get_intensity_percentile(0.75, open(wig, 'r'))
    annotations = ["Promoter", "Enhancer", "Weak"]
    length = len(annotations)
    ground_truth_list = []
    predicted_list = []
      
    #Get precision and recall for each type.
    [pred, gt] = get_labels_and_ground_truth(bed, sig, wig, annotations, threshold)
    precision = dict()
    recall = dict()
    for i in range(length):
        precision[i] = precision_score(gt[:, i], pred[:, i])
        recall[i] = recall_score(gt[:, i], pred[:, i])
    #Get precision and recall for TSS-based promoters.
    #[pred_promoter, gt_promoter] = get_tss_labels_and_ground_truth(tss_bed, tss_sig, wig, ["Promoter", "Not_Promoter"], threshold) 
    #precision["tss"] = precision_score(gt_promoter[:, 0], pred_promoter[:, 0])
    #recall["tss"] = recall_score(gt_promoter[:, 0], pred_promoter[:, 0])
    
    #Get precision and recall for combined annotations (OR).
    # [pred_or, gt_or] = get_labels_and_ground_truth(or_bed, or_sig, wig, annotations, threshold)
    # precision_or = dict()
    # recall_or = dict()
    # for i in range(length):
        # precision_or[i] = precision_score(gt_or[:, i], pred_or[:, i])
        # recall_or[i] = recall_score(gt_or[:, i], pred_or[:, i])
        
    #Get precision and recall for combined annotations (AND).
    [pred_and, gt_and] = get_labels_and_ground_truth(and_bed, and_sig, wig, annotations, threshold)
    precision_and = dict()
    recall_and = dict()
    for i in range(length):
        precision_and[i] = precision_score(gt_and[:, i], pred_and[:, i])
        recall_and[i] = recall_score(gt_and[:, i], pred_and[:, i])
        
    # #Get precision and recall for permuted annotations.
    # # [pred_perm, gt_perm] = get_labels_and_ground_truth(perm_bed, perm_sig, wig, annotations, threshold)
    # # precision_perm = dict()
    # # recall_perm = dict()
    # # for i in range(length):
        # # precision_perm[i] = precision_score(gt_perm[:, i], pred_perm[:, i])
        # # recall_perm[i] = recall_score(gt_perm[:, i], pred_perm[:, i])
        
    # #Get precision and recall for RPKM annotations.
    # [pred_rpkm, gt_rpkm] = get_labels_and_ground_truth(rpkm_bed, rpkm_sig, wig, annotations, threshold)
    # precision_rpkm = dict()
    # recall_rpkm = dict()
    # for i in range(length):
        # precision_rpkm[i] = precision_score(gt_rpkm[:, i], pred_rpkm[:, i])
        # recall_rpkm[i] = recall_score(gt_rpkm[:, i], pred_rpkm[:, i])
        
    #return [precision, recall, precision_or, recall_or, precision_perm, recall_perm, precision_rpkm, recall_rpkm, pred.shape[0], threshold, pred, gt]
    return [precision, recall, precision_and, recall_and, pred.shape[0], threshold, pred, gt]
    
"""
Print precision and recall for all chromosomes.
""" 
#def print_report(precision, recall, precision_or, recall_or, precision_and, recall_and, precision_perm, recall_perm, precision_rpkm, recall_rpkm, chrom, win, cell, pr):
def print_report(precision, recall, precision_and, recall_and, chrom, win, cell, pr):
    #Print the cell line, chromosome, and window information
    #Weak   (Shape, TSS + Shape, Perm)
    report = open(pr + "/report_" + chrom, "w")
    report.write(cell + "\n")
    report.write(chrom + "\n")
    report.write(win + "\n\n")
    
    #All shape-based predictions
    for i in range(0, 3):
        report.write(str(precision[i]) + "\t" + str(recall[i]) + "\n")
    report.write("\n")
    
    # #TSS predictions for promoters
    # report.write(str(precision["tss"]) + "\t" + str(recall["tss"]) + "\n")
    # report.write("\n")
    
    # #TSS OR Shape predictions
    # for i in range(0, 3):
        # report.write(str(precision_or[i]) + "\t" + str(recall_or[i]) + "\n")
    # report.write("\n")
    
    #TSS AND Shape predictions
    for i in range(0, 3):
        report.write(str(precision_and[i]) + "\t" + str(recall_and[i]) + "\n")
    report.write("\n")
    
    # #Permuted predictions
    # for i in range(0, 3):
        # report.write(str(precision_perm[i]) + "\t" + str(recall_perm[i]) + "\n")
        
    # #Permuted predictions
    # for i in range(0, 3):
        # report.write(str(precision_rpkm[i]) + "\t" + str(recall_rpkm[i]) + "\n")
    
    #Close the report.
    report.close()

#Get the percentage of the chromosome belonging to each ChromHMM annotation.
#def save_scatterplot(our_precision, our_recall, tss_precision, tss_recall, or_precision, or_recall, perm_precision, perm_recall, rpkm_precision, rpkm_recall, out, cell, src, indices_to_highlight, chroms, separate_legend):
#def save_scatterplot(our_precision, our_recall, tss_precision, tss_recall, or_precision, or_recall, and_precision, and_recall, rpkm_precision, rpkm_recall, out, cell, src, indices_to_highlight, chroms, separate_legend):
def save_scatterplot(our_precision, our_recall, and_precision, and_recall, out, cell, src, indices_to_highlight, chroms, separate_legend):

    #Set colors and symbols for plotting.
    enhancer_color = "blue"
    promoter_color = "green"
    weak_color = "red"
    our_symbol = "*"
    tss_symbol = "D"
    or_symbol = "X"
    and_symbol = "P"
    perm_symbol = "s"
    rpkm_symbol = "|"
    our_size = 10
    tss_size = 5
    plus_size = 7
    perm_size = 5
    rpkm_size = 5
    
    #Set the axes, title, and maximum.
    plt.ylim(-0.05,1.05)
    plt.xlim(-0.05,1.05)
    plt.title("Precision and Recall for " + cell + " Predicted from " + src)
    plt.xlabel("Precision")
    plt.ylabel("Recall")
    
    #Plot our data.
    plt.scatter(our_precision[0,:], our_recall[0,:], c = promoter_color, marker = our_symbol, edgecolor = "black", s = 70)
    plt.scatter(our_precision[1,:], our_recall[1,:], c = enhancer_color, marker = our_symbol, edgecolor = "black", s = 70)
    plt.scatter(our_precision[2,:], our_recall[2,:], c = weak_color, marker = our_symbol, edgecolor = "black", s = 70)
    
    #Add text for chromosomes we learned on.
    for i in indices_to_highlight:
        plt.text(our_precision[0,i] + 0.01, our_recall[0,i] + 0.01, chroms[i], size = 8)
        plt.text(our_precision[1,i] + 0.01, our_recall[1,i] + 0.01, chroms[i], size = 8)
        plt.text(our_precision[2,i] + 0.01, our_recall[2,i] + 0.01, chroms[i], size = 8)
        
    # #Plot TSS data.
    # plt.scatter(tss_precision[0,:], tss_recall[0,:], c = promoter_color, marker = tss_symbol)
    
    # #Plot combined (OR) data.
    # plt.scatter(or_precision[0,:], or_recall[0,:], c = promoter_color, marker = or_symbol)
    # plt.scatter(or_precision[1,:], or_recall[1,:], c = enhancer_color, marker = or_symbol)
    # plt.scatter(or_precision[2,:], or_recall[2,:], c = weak_color, marker = or_symbol)
    
    #Plot combined (AND) data.
    plt.scatter(and_precision[0,:], and_recall[0,:], c = promoter_color, marker = and_symbol)
    plt.scatter(and_precision[1,:], and_recall[1,:], c = enhancer_color, marker = and_symbol)
    plt.scatter(and_precision[2,:], and_recall[2,:], c = weak_color, marker = and_symbol)
    
    # #Plot permuted data.
    # #plt.scatter(perm_precision[0,:], perm_recall[0,:], c = "white", marker = perm_symbol, edgecolor = promoter_color)
    # #plt.scatter(perm_precision[1,:], perm_recall[1,:], c = "white", marker = perm_symbol, edgecolor = enhancer_color)
    # #plt.scatter(perm_precision[2,:], perm_recall[2,:], c = "white", marker = perm_symbol, edgecolor = weak_color)
    
    # #Plot rpkm data.
    # plt.scatter(rpkm_precision[0,:], rpkm_recall[0,:], c = promoter_color, marker = rpkm_symbol)
    # plt.scatter(rpkm_precision[1,:], rpkm_recall[1,:], c = enhancer_color, marker = rpkm_symbol)
    # plt.scatter(rpkm_precision[2,:], rpkm_recall[2,:], c = weak_color, marker = rpkm_symbol)
    
    if separate_legend:
        plt.savefig(out + "precision_or_recall_nolegend" + src + ".png")
        plt.close()
    
    #Add the legend.
    legend_elements = [Patch(facecolor=enhancer_color, label='Enhancer'),
                        Patch(facecolor=promoter_color, label='Promoter'),
                        Patch(facecolor=weak_color, label='Weak'),
                        Line2D([0], [0], marker=our_symbol, markerfacecolor='black', label='Shape-Based', color = 'white',
                           markersize=our_size),
                        Line2D([0], [0], marker=tss_symbol, markerfacecolor='black', label='TSS-Based', color = 'white',
                           markersize=tss_size),
                        Line2D([0], [0], marker=or_symbol, markerfacecolor='black', label='Shape OR TSS', color = 'white',
                           markersize=plus_size),
                        Line2D([0], [0], marker=and_symbol, markerfacecolor='black', label='Shape AND TSS', color = 'white',
                           markersize=plus_size),
                        #Line2D([0], [0], marker=perm_symbol, markerfacecolor="None", markeredgecolor = 'black', color = 'white', label='Permuted', markersize=perm_size),
                        Line2D([0], [0], marker=rpkm_symbol, markerfacecolor='black', markeredgecolor = 'black', color = 'white', label='RPKM-Based',
                           markersize=rpkm_size),
                       ]
    #fig, ax = plt.subplots()
    plt.legend(handles=legend_elements, loc="lower right")
    
    #Save the plot.
    plt.savefig(out + "precision_or_recall.png")
    plt.close()
    
#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def get_labels_and_ground_truth(bed_file, sig_file, wig, annotations, threshold):

    #Set up percentage matrix.
    vec_pred = list()
    vec_gt = list()
    
    #Get scores and labels for each bed file.
    bed = np.genfromtxt(bed_file , delimiter='\t', dtype = str)
    sigf = open(sig_file, "r")

    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    prev_start = -1
    prev_end = -1
    sig_i = 0
    sig = [float(s) for s in sigf.readline().split(",")]
    sum_vec = np.zeros(3)
    
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
            if sig_i != 0:
                sig_s = sigf.readline()
                sig = [float(s) for s in sig_s.split(",")]
            sum_vec = np.zeros(3)
            sig_i += 1
       
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
        # elif a == "13_ReprPC" or a == "ReprPCWk":
            # if total_peak_size > 0:
                # sum_vec[3] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            # else:
                # sum_vec[3] += anno_length
        elif a == "9_Het" or a == "15_Quies":
            if total_peak_size > 0:
                sum_vec[2] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
            else:
                sum_vec[2] += anno_length
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
            sum_vec = np.zeros(4)
            
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
            # elif a == "13_ReprPC" or a == "ReprPCWk":
                # if total_peak_size > 0:
                    # sum_vec[3] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                # else:
                    # sum_vec[3] += anno_length
            elif a == "9_Het" or a == "15_Quies":
                if total_peak_size > 0:
                    sum_vec[2] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[2] += anno_length
            elif anno_length != 0:
                if total_peak_size > 0:
                    sum_vec[3] += count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end)
                else:
                    sum_vec[3] += anno_length
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
                if sum_vec[0] == 0 and sum_vec[1] == 0 and sum_vec[2] == 0 and sum_vec[3] == 0:
                    del vec_gt[len(vec_gt) - 1]
                #Set count and unannotated count to 0.
                not_annotated_count = 0
                count_in_region = 0

        #Get the previous data, if applicable.
        prev_start = current_start
        prev_end = current_start
            
    #Write a report with the percentage of unknowns and, of the unknowns, the percentage of each annotation type.
    out = open(out_path + "unknown_dist_" + chrom, "w")
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
    percent_other = count_all[3] / ground_truth.shape[0]
    out.write(str(percent_other) + "\n")
    
#Plot a heatmap of the mispredictions.
def save_misprediction_heatmap(predictions, ground_truth, path, cell, src):
    
    #Annotation names
    names = ["Promoter", "Enhancer", "Weak"]

    #Find the counts of each misprediction.
    where_promoters = np.where(ground_truth[:,0] == 1)[0]
    promoters = predictions[where_promoters]
    promoter_pred_vec = np.sum(promoters, axis = 0) / promoters.shape[0]
    where_enhancers = np.where(ground_truth[:,1] == 1)[0]
    enhancers = predictions[where_enhancers]
    enhancer_pred_vec = np.sum(enhancers, axis = 0) / enhancers.shape[0]
    where_weaks = np.where(ground_truth[:,2] == 1)[0]
    weaks = predictions[where_weaks]
    weak_pred_vec = np.sum(weaks, axis = 0) / weaks.shape[0]
    mispredictions = np.stack([promoter_pred_vec, enhancer_pred_vec, weak_pred_vec])
    
    #Plot the heatmap.
    f = plt.figure(1, figsize=(8, 5))
    heatmap = sns.heatmap(mispredictions, cbar=True, cmap="BuGn", fmt="d", vmin = 0, vmax = 1, xticklabels = names, yticklabels = names)
    plt.title("Composition of Predictions per ChromHMM Regulatory Category from " + src + " to " + cell)
    plt.xlabel("Predicted Regulatory Category")
    plt.ylabel("ChromHMM Regulatory Category")
    plt.yticks(rotation = 0)
    fig = heatmap.get_figure()
    fig.savefig(path + "misprediction_heatmap" + src + ".png")
 
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