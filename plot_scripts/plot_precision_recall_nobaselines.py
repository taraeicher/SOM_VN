import numpy as np
import scipy as sp
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import glob, os
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
matplotlib.rcParams.update({'font.size': 24})

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the chromHMM and annotated files.
    plot_out = sys.argv[3]
    title = sys.argv[4]
    all_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
    precision_all = np.zeros((3, len(all_chroms)))
    recall_all = np.zeros((3, len(all_chroms)))
    predictions_all = []
    ground_truth_all = []

    #Plot ROC curve for each chromosome. Plot curve for each annotation separately.
    c = 0
    for chrom in all_chroms:
        #Get files for each category.
        our_bed = sys.argv[1] + "anno" + str(chrom) + ".bed"
        our_sig = sys.argv[2] + "clusters_anno" + str(chrom)
        wig = sys.argv[5] + str(chrom) + ".wig"
        
        #Get all precision and recall values.
        [precision, recall, total, threshold, predictions, ground_truth, fpr] = get_all_precision_and_recall(our_bed, our_sig, wig, chrom)
        predictions_all.append(predictions)
        ground_truth_all.append(np.asarray(ground_truth))
        
        for i in range(0,3):
            precision_all[i,c] = precision[i]
            recall_all[i,c] = recall[i]

        c += 1
        
        #Save a scatterplot.
        save_scatterplot(precision_all, recall_all, plot_out, title, all_chroms)
    
#Plot the ROC curve based on ground truth and prediction for each cutoff. Plot separate lines for each
#annotation type. Consolidate all chromosomes.
def get_all_precision_and_recall(bed, sig, wig, chrom):

    #Get actual annotation and ground truth for all annotations and for all unannotated regions.
    threshold = wsu.get_intensity_percentile(0.75, open(wig, 'r'), 0)
    annotations = ["Promoter", "Enhancer", "Weak"]
    length = len(annotations)
    ground_truth_list = []
    predicted_list = []
      
    #Get precision and recall for each type.
    [pred, gt] = get_labels_and_ground_truth(bed, sig, wig, annotations, threshold)
    precision = dict()
    recall = dict()
    fpr = dict()
    for i in range(length):
        precision[i] = precision_score(gt[:, i], pred[:, i])
        recall[i] = recall_score(gt[:, i], pred[:, i])
        fp = len(np.where((pred[:, i] == 1) & (gt[:, i] == 0))[0])
        tn = len(np.where((pred[:, i] == 0) & (gt[:, i] == 0))[0])
        fpr[i] = fp / (fp + tn)
    return [precision, recall, pred.shape[0], threshold, pred, gt, fpr]
    

#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def save_scatterplot(our_precision, our_recall, out, title, indices_to_highlight):

    #Set colors and symbols for plotting.
    enhancer_color = "gray"
    promoter_color = "black"
    weak_color = "white"
    weak_outline_color = "silver"
    edge_weak = "black"
    our_symbol = "*"
    our_size = 10
    factor = 40
    plt.gcf().subplots_adjust(bottom=0.20)
    plt.gcf().subplots_adjust(left=0.20)
    
    #Set the axes, title, and maximum.
    plt.ylim(-0.05,1.05)
    plt.xlim(-0.05,1.05)
    plt.title(title)
    plt.xlabel("Precision")
    plt.ylabel("Recall")
    
    #Plot our data.
    plt.scatter(our_precision[0,:], our_recall[0,:], c = promoter_color, marker = our_symbol, edgecolor = "black", s = our_size * factor)
    plt.scatter(our_precision[1,:], our_recall[1,:], c = enhancer_color, marker = our_symbol, edgecolor = "black", s = our_size * factor)
    plt.scatter(our_precision[2,:], our_recall[2,:], c = weak_color, marker = our_symbol, edgecolor = "black", s = our_size * factor)

    #Save
    plt.savefig(out + "precision_recall_nolegend" + title + ".png")
    plt.close()
    
    #Add the legend.
    legend_elements = [Patch(facecolor=enhancer_color, label='Enhancer'),
                        Patch(facecolor=promoter_color, label='Promoter'),
                        Patch(facecolor=weak_color, edgecolor = edge_weak, label='Weak')
                       ]
    plt.legend(handles=legend_elements, loc="lower right")
    
    #Save the plot.
    plt.savefig(out + "legend.png")
    plt.close()
    
#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def get_labels_and_ground_truth(bed_file, sig_file, wig, annotations, threshold):
    
    #Set up percentage matrix.
    vec_pred = list()
    vec_gt = list()
    final_stack_pred = np.empty((0, 0))
    final_stack_gt = np.empty((0, 0))

    #Get scores and labels for each bed file.
    bed = np.genfromtxt(bed_file , delimiter='\t', dtype = str)
    sigf = open(sig_file, "r")
    
    #Loop through bed file to compute percentage for each region.
    current_start = -1
    current_end = -1
    prev_start = -1
    prev_end = -1
    sig_i = 0
    #Do not move forward if the first line is blank in the sig file.
    try:
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
            total_peak_size = wsu.count_above(threshold, "", sig, current_start, current_end, current_start, current_end, BIN_SIZE)
            if a == "1_TssA" or a == "2_TssAFlnk" or a == "10_TssBiv" or a == "11_BivFlnk":
                if total_peak_size > 0:
                    sum_vec[0] += wsu.count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end, BIN_SIZE)
                else:
                    sum_vec[0] += anno_length
            elif a == "6_EnhG" or a == "7_Enh" or a == "12_EnhBiv":
                if total_peak_size > 0:
                    sum_vec[1] += wsu.count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end, BIN_SIZE)
                else:
                    sum_vec[1] += anno_length
            elif a == "9_Het" or a == "15_Quies":
                if total_peak_size > 0:
                    sum_vec[2] += wsu.count_above(threshold, a, sig, current_start, current_end, anno_start, anno_end, BIN_SIZE)
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
            
        #Stack all values.
        final_stack_pred = np.stack(vec_pred)
        final_stack_gt = np.stack(vec_gt)
    except:
        pass
    #Return value.
    return [final_stack_pred, final_stack_gt]
    
if __name__ == "__main__":
    main()