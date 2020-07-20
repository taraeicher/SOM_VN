import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import glob, os
import math
matplotlib.rcParams.update({'font.size': 16})

"""
This program takes the precision and recall values from all chromosomes and all runs and plots them for weak, promoter, and enhancer in a hexbin plot.
"""
def main():

    #Get inputs.
    cell = ["A549", "Brain", "H1"]
    iters = int(sys.argv[1])
    filename = sys.argv[2]
    
    #Read all precision and recall and save to a hexbin plot.
    num_plots = 3
    for i in range(num_plots):
        [precision_weak, recall_weak, precision_promoter, recall_promoter, precision_enhancer, recall_enhancer] = read_all_precision_and_recall(cell[i], iters)
        print("Mean precision weak: " + cell[i] + " " + str(np.mean(precision_weak)))
        print("Mean recall weak: " + cell[i] + " " + str(np.mean(recall_weak)))
        print("Mean precision promoter: " + cell[i] + " " + str(np.mean(precision_promoter)))
        print("Mean recall promoter: " + cell[i] + " " + str(np.mean(recall_promoter)))
        print("Mean precision enhancer: " + cell[i] + " " + str(np.mean(precision_enhancer)))
        print("Mean recall enhancer: " + cell[i] + " " + str(np.mean(recall_enhancer)))

        plt.subplot(2,2,i+1)
        make_hexbin_plots(i,precision_weak, recall_weak, precision_promoter, recall_promoter, precision_enhancer, recall_enhancer, cell[i])
    
    #Save the plots.
    plt.savefig(filename + ".png", bbox_inches = "tight", dpi = 300)

    legend_elements = [Patch(facecolor="green", label='Promoter'),
                        Patch(facecolor="blue", label='Enhancer'),
                        Patch(facecolor="red", label='Weak')
                       ]
    plt.close()
    plt.legend(handles=legend_elements, loc="lower right")
    plt.axis("off")
    plt.savefig(filename + "_legend.png", bbox_inches = "tight", dpi = 300)

"""
Read precision and recall values from each regulatory category from each chromosome from each run.
"""
def read_all_precision_and_recall(cell, iter_count):
    
    #Declare the lists to fill in.
    pw = []
    rw = []
    pp = []
    rp = []
    pe = []
    re = []
    
    #For each iteration, fill in the values.
    for i in range(1, iter_count + 1):
        directory = "/fs/scratch/PAS0272/DNase_SOM_06042020/DNase_SOM/old_results/" + cell + "/precision_recall_" + str(i)
        
        #For each chromosome in the directory for an iteration, fill in the values.
        for filename in os.listdir(directory):
            f = open(os.path.join(directory, filename), "r") 
            
            #Skip first four lines.
            for j in range(0, 4):
                junk = f.readline()
                
            #Structure: precision\trecall.
            #Order: Promoter, Enhancer, Weak
            pr_pair = f.readline().split("\t")
            pp.append(float(pr_pair[0]))
            rp.append(float(pr_pair[1]))
            pr_pair = f.readline().split("\t")
            pe.append(float(pr_pair[0]))
            re.append(float(pr_pair[1]))
            pr_pair = f.readline().split("\t")
            pw.append(float(pr_pair[0]))
            rw.append(float(pr_pair[1]))
            
    #Return lists.
    return [pw, rw, pp, rp, pe, re]
                      
"""
Save overlapping hexbin plots of precision and recall for all regulatory categories.
"""
def make_hexbin_plots(pos,pw, rw, pp, rp, pe, re, cell):
    
    #Make plots for each regulatory category.
    sns.kdeplot(rw, pw, cmap=plt.cm.Reds)
    nonzeros = (np.nonzero(pe)[0]).tolist()
    zeros = (np.where(np.asarray(pe) == 0.0))
    try:
        sns.kdeplot(re, pe, cmap=plt.cm.Blues)
    except:
        sns.kdeplot(np.asarray(re)[nonzeros], np.asarray(pe)[nonzeros], cmap=plt.cm.Blues)
        sns.regplot(x=np.array([0]), y=np.array([0]), scatter=True, fit_reg=False, marker='o',
            scatter_kws={"s": 10}, color = "black") 
    sns.kdeplot(pp, rp, cmap=plt.cm.Greens)
    
    #Set x and y axis names.
    if pos > 0:
       plt.xlabel("Recall")
    if pos != 1:
    	plt.ylabel("Precision")
    if pos == 0:
        plt.xticks([])
    if pos == 1:
        plt.yticks([])
    plt.title(cell)
    
    
if __name__ == "__main__":
    main()
