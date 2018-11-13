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

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Set up the files and source-dest combinations.
    files = []
    combo = []
    cells = ["A549", "Brain", "H1"]
    for src in cells:
        for dest in cells:
            files.append(sys.argv[1] + dest + sys.argv[2] + src)
            combo.append(src[0] + "-" + dest[0])
    plot_out = sys.argv[4]
    all_chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21']
       
    #Set up precision and recall containers.
    precision_all = np.zeros((3, 9, len(all_chroms)))
    recall_all = np.zeros((3, 9, len(all_chroms)))
    precision_tss_all = np.zeros((9, len(all_chroms)))
    recall_tss_all = np.zeros((9, len(all_chroms)))
    precision_distal_all = np.zeros((9, len(all_chroms)))
    recall_distal_all = np.zeros((9, len(all_chroms)))
    precision_or_all = np.zeros((3, 9, len(all_chroms)))
    recall_or_all = np.zeros((3, 9, len(all_chroms)))
    precision_and_all = np.zeros((3, 9, len(all_chroms)))
    recall_and_all = np.zeros((3, 9, len(all_chroms)))
    precision_perm_all = np.zeros((3, 9, len(all_chroms)))
    recall_perm_all = np.zeros((3, 9, len(all_chroms)))
    precision_rpkm_all = np.zeros((3, 9, len(all_chroms)))
    recall_rpkm_all = np.zeros((3, 9, len(all_chroms)))

    #Plot ROC curve for each chromosome. Plot curve for each annotation separately.
    f = 0
    for file in files:
        c = 0
        for chrom in all_chroms:
            #Get files for each category.
            this_file = open(file + sys.argv[3] + chrom)
            
            #Get all precision and recall values.
            [precision, recall, precision_tss, recall_tss, precision_distal, recall_distal, precision_or, recall_or, precision_and, recall_and, precision_perm, recall_perm, precision_rpkm, recall_rpkm] = read_pr_from_file(this_file)
            
            for i in range(0,3):
                precision_all[i,f,c] = precision[i]
                recall_all[i,f,c] = recall[i]
                precision_or_all[i,f,c] = precision_or[i]
                recall_or_all[i,f,c] = recall_or[i]
                precision_and_all[i,f,c] = precision_and[i]
                recall_and_all[i,f,c] = recall_and[i]
                precision_perm_all[i,f,c] = precision_perm[i]
                recall_perm_all[i,f,c] = recall_perm[i]
                precision_rpkm_all[i,f,c] = precision_rpkm[i]
                recall_rpkm_all[i,f,c] = recall_rpkm[i]
                
            precision_tss_all[f,c] = precision_tss
            recall_tss_all[f,c] = recall_tss
            precision_distal_all[f,c] = precision_distal
            recall_distal_all[f,c] = recall_distal
            c += 1
        f += 1
        
    #Sum over chromosomes.
    precision = np.sum(precision_all, axis = 2) / len(all_chroms)
    recall = np.sum(recall_all, axis = 2) / len(all_chroms)
    p_tss = np.sum(precision_tss_all, axis = 1) / len(all_chroms)
    r_tss = np.sum(recall_tss_all, axis = 1) / len(all_chroms)
    p_distal = np.sum(precision_distal_all, axis = 1) / len(all_chroms)
    r_distal = np.sum(recall_distal_all, axis = 1) / len(all_chroms)
    precision_or = np.sum(precision_or_all, axis = 2) / len(all_chroms)
    recall_or = np.sum(recall_or_all, axis = 2) / len(all_chroms)
    precision_and = np.sum(precision_and_all, axis = 2) / len(all_chroms)
    recall_and = np.sum(recall_and_all, axis = 2) / len(all_chroms)
    precision_perm = np.sum(precision_perm_all, axis = 2) / len(all_chroms)
    recall_perm = np.sum(recall_perm_all, axis = 2) / len(all_chroms)
    precision_rpkm = np.sum(precision_rpkm_all, axis = 2) / len(all_chroms)
    recall_rpkm = np.sum(recall_rpkm_all, axis = 2) / len(all_chroms)
    
    #Save a scatterplot with all precision and recall values.
    save_scatterplot(precision, recall, p_tss, r_tss, p_distal, r_distal, precision_or, recall_or, precision_and, recall_and, precision_perm, recall_perm, precision_rpkm, recall_rpkm, plot_out, combo)
    
"""
Return the precision and recall values from a file.
""" 
def read_pr_from_file(report):
    #Open and read the report.
    junk = report.readline()
    junk = report.readline()
    junk = report.readline()
    junk = report.readline()
    
    #All shape-based predictions
    precision = np.zeros(3)
    recall = np.zeros(3)
    for i in range(0, 3):
        line = report.readline().split("\t")
        precision[i] = line[0]
        recall[i] = line[1]
        
    #TSS predictions for promoters
    junk = report.readline()
    line = report.readline().split("\t")
    precision_tss = line[0]
    recall_tss = line[1]
    
    #TSS predictions for distal
    junk = report.readline()
    line = report.readline().split("\t")
    precision_distal = line[0]
    recall_distal = line[1]
    
    #TSS OR Shape predictions
    precision_or = np.zeros(3)
    recall_or = np.zeros(3)
    junk = report.readline()
    for i in range(0, 3):
        line = report.readline().split("\t")
        precision_or[i] = line[0]
        recall_or[i] = line[1]
        
    #TSS AND Shape predictions
    precision_and = np.zeros(3)
    recall_and = np.zeros(3)
    junk = report.readline()
    for i in range(0, 3):
        line = report.readline().split("\t")
        precision_and[i] = line[0]
        recall_and[i] = line[1]
    
    #Permuted predictions
    precision_perm = np.zeros(3)
    recall_perm = np.zeros(3)
    #junk = report.readline()
    # for i in range(0, 3):
        # line = report.readline().split("\t")
        # precision_perm[i] = line[0]
        # recall_perm[i] = line[1]
        
    #RPKM predictions
    precision_rpkm = np.zeros(3)
    recall_rpkm = np.zeros(3)
    junk = report.readline()
    for i in range(0, 3):
        line = report.readline().split("\t")
        precision_rpkm[i] = line[0]
        recall_rpkm[i] = line[1]
    
    #Return all values.
    return [precision, recall, precision_tss, recall_tss, precision_distal, recall_distal, precision_or, recall_or, precision_and, recall_and, precision_perm, recall_perm, precision_rpkm, recall_rpkm]

#Get the percentage of the chromosome belonging to each ChromHMM annotation.
def save_scatterplot(our_precision, our_recall, tss_precision, tss_recall, distal_precision, distal_recall, or_precision, or_recall, and_precision, and_recall, perm_precision, perm_recall, rpkm_precision, rpkm_recall, out, labels):

    #Set colors and symbols for plotting.
    enhancer_color = "gray"
    promoter_color = "black"
    distal_color = "darkviolet"
    weak_color = "silver"
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
    plt.xlabel("Precision")
    plt.ylabel("Recall")
    
    #Plot our data.
    plt.scatter(our_precision[0,:], our_recall[0,:], c = promoter_color, marker = our_symbol, edgecolor = "black", s = 70)
    enhancer_precision = [x for i,x in enumerate(our_precision[1,:]) if (labels[i] != "H-H" and labels[i] != "H-B" and labels[i] != "H-A")]
    enhancer_recall= [x for i,x in enumerate(our_recall[1,:]) if (labels[i] != "H-H" and labels[i] != "H-B" and labels[i] != "H-A")]
    plt.scatter(enhancer_precision, enhancer_recall, c = enhancer_color, marker = our_symbol, edgecolor = "black", s = 70)
    plt.scatter(our_precision[2,:], our_recall[2,:], c = weak_color, marker = our_symbol, edgecolor = "black", s = 70)
    
    #Plot AND data.
    plt.scatter(and_precision[0,:], and_recall[0,:], c = promoter_color, marker = and_symbol)
    and_enhancer_precision = [x for i,x in enumerate(and_precision[1,:]) if (labels[i] != "H-H" and labels[i] != "H-B" and labels[i] != "H-A")]
    and_enhancer_recall= [x for i,x in enumerate(and_recall[1,:]) if (labels[i] != "H-H" and labels[i] != "H-B" and labels[i] != "H-A")]
    plt.scatter(and_enhancer_precision, and_enhancer_recall, c = enhancer_color, marker = and_symbol)
    plt.scatter(and_precision[2,:], and_recall[2,:], c = weak_color, marker = and_symbol)
    
    #Add labels for our experiments.
    for i in range(0, len(labels)):
        plt.text(our_precision[0,i] + 0.01, our_recall[0,i] + 0.01, labels[i], size = 8)
        if labels[i] != "H-H" and labels[i] != "H-B" and labels[i] != "H-A" and labels[i] != "B-A":
            plt.text(our_precision[1,i] + 0.01, our_recall[1,i] + 0.01, labels[i], size = 8)
        elif labels[i] == "B-A":
            plt.text(our_precision[1,i] + 0.02, our_recall[1,i] - 0.01, labels[i], size = 8)
        if labels[i] != "B-H":
            plt.text(our_precision[2,i] + 0.01, our_recall[2,i] + 0.01, labels[i], size = 8)
        else:
            plt.text(our_precision[2,i] + 0.02, our_recall[2,i], labels[i], size = 8)
        
    #Plot TSS data.
    plt.scatter(tss_precision, tss_recall, c = promoter_color, marker = tss_symbol)
    #plt.scatter(distal_precision, distal_recall, c = distal_color, marker = tss_symbol)
    
    #Add labels for our experiments.
    plt.text(tss_precision[0] + 0.01, tss_recall[0] + 0.01, "A549", size = 8)
    plt.text(tss_precision[1] + 0.01, tss_recall[1] + 0.01, "Brain", size = 8)
    plt.text(tss_precision[2] + 0.02, tss_recall[2], "H1", size = 8)
    
    #Plot combined data.
    # plt.scatter(or_precision[0,:], or_recall[0,:], c = promoter_color, marker = or_symbol)
    # plt.scatter(or_precision[1,:], or_recall[1,:], c = enhancer_color, marker = or_symbol)
    # plt.scatter(or_precision[2,:], or_recall[2,:], c = weak_color, marker = or_symbol)
    
    # #Add labels for our experiments.
    # for i in range(0, len(labels)):
        # plt.text(or_precision[0,i] + 0.01, or_recall[0,i] + 0.01, labels[i], size = 8)
        # plt.text(or_precision[1,i] + 0.01, or_recall[1,i] + 0.01, labels[i], size = 8)
        # plt.text(or_precision[2,i] + 0.01, or_recall[2,i] + 0.01, labels[i], size = 8)
    
    #Plot combined data.
    # plt.scatter(and_precision[0,:], and_recall[0,:], c = promoter_color, marker = or_symbol)
    # plt.scatter(and_precision[1,:], and_recall[1,:], c = enhancer_color, marker = or_symbol)
    # plt.scatter(and_precision[2,:], and_recall[2,:], c = weak_color, marker = or_symbol)
    
    # #Add labels for our experiments.
    # for i in range(0, len(labels)):
        # plt.text(and_precision[0,i] + 0.01, and_recall[0,i] + 0.01, labels[i], size = 8)
        # plt.text(and_precision[1,i] + 0.01, and_recall[1,i] + 0.01, labels[i], size = 8)
        # plt.text(and_precision[2,i] + 0.01, and_recall[2,i] + 0.01, labels[i], size = 8)
        
    #Plot permuted data.
    # plt.scatter(perm_precision[0,:], perm_recall[0,:], c = "white", marker = perm_symbol, edgecolor = promoter_color)
    # plt.scatter(perm_precision[1,:], perm_recall[1,:], c = "white", marker = perm_symbol, edgecolor = enhancer_color)
    # plt.scatter(perm_precision[2,:], perm_recall[2,:], c = "white", marker = perm_symbol, edgecolor = weak_color)
    
    # #Add labels for our experiments.
    # for i in range(0, len(labels)):
        # plt.text(perm_precision[0,i] + 0.01, perm_recall[0,i] + 0.01, labels[i], size = 8)
        # plt.text(perm_precision[1,i] + 0.01, perm_recall[1,i] + 0.01, labels[i], size = 8)
        # plt.text(perm_precision[2,i] + 0.01, perm_recall[2,i] + 0.01, labels[i], size = 8)
    
    #Plot RPKM data.
    plt.scatter(rpkm_precision[0,:], rpkm_recall[0,:], c = promoter_color, marker = rpkm_symbol)
    plt.scatter(rpkm_precision[1,:], rpkm_recall[1,:], c = enhancer_color, marker = rpkm_symbol)
    plt.scatter(rpkm_precision[2,:], rpkm_recall[2,:], c = weak_color, marker = rpkm_symbol)
    
    #Add labels for our experiments.
    plt.text(rpkm_precision[0,0] + 0.01, rpkm_recall[0,0], "A549", size = 8)
    plt.text(rpkm_precision[0,1] + 0.01, rpkm_recall[0,1] + 0.01, "Brain", size = 8)
    plt.text(rpkm_precision[0,2] + 0.01, rpkm_recall[0,2] - 0.01, "H1", size = 8)
    plt.text(rpkm_precision[1,0] + 0.01, rpkm_recall[1,0] + 0.01, "A549", size = 8)
    plt.text(rpkm_precision[1,1] + 0.01, rpkm_recall[1,1] + 0.01, "Brain", size = 8)
    plt.text(rpkm_precision[1,2] + 0.01, rpkm_recall[1,2] + 0.01, "H1", size = 8)
    plt.text(rpkm_precision[2,0] + 0.01, rpkm_recall[2,0] + 0.01, "A549", size = 8)
    plt.text(rpkm_precision[2,1] + 0.01, rpkm_recall[2,1] + 0.01, "Brain", size = 8)
    plt.text(rpkm_precision[2,2] + 0.01, rpkm_recall[2,2] + 0.01, "H1", size = 8)
    
    #Close and put legend in separate file.
    plt.savefig(out + "precision_and_recall_all.png")
    plt.close()
    
    #Add the legend.
    legend_elements = [Patch(facecolor=enhancer_color, label='Enhancer'),
                        Patch(facecolor=promoter_color, label='Promoter'),
                        Patch(facecolor=weak_color, label='Weak'),
                        #Patch(facecolor=distal_color, label='TSS Distal'),
                        Line2D([0], [0], marker=our_symbol, markerfacecolor='black', color = 'white', label='Shape-Based',
                           markersize=our_size),
                        Line2D([0], [0], marker=tss_symbol, markerfacecolor='black', color = 'white', label='TSS-Based',
                           markersize=tss_size),
                        #Line2D([0], [0], marker=or_symbol, markerfacecolor='black', color = 'white', label='Shape + TSS',
                           #markersize=plus_size),
                        Line2D([0], [0], marker=and_symbol, markerfacecolor='black', color = 'white', label='Shape + TSS',
                           markersize=plus_size),
                        #Line2D([0], [0], marker=perm_symbol, markerfacecolor="None", markeredgecolor = 'black', color = 'white', label='Permuted', markersize=perm_size),
                        Line2D([0], [0], marker=rpkm_symbol, markerfacecolor='black', markeredgecolor = 'black', color = 'white', label='RPKM-Based', markersize=rpkm_size),
                       ]
    #fig, ax = plt.subplots()
    plt.legend(handles=legend_elements, loc="lower right")
    
    #Save the plot.
    plt.savefig(out + "legend.png")
    plt.close()            
 
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