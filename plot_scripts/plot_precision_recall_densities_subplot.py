import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import glob, os
matplotlib.rcParams.update({'font.size': 12})

"""
This program takes the precision and recall values from all chromosomes and all runs and plots them for weak, promoter, and enhancer in a hexbin plot.
"""
def main():

    #Get inputs.
    cell = ["A549", "Brain", "H1"]
    iters = int(sys.argv[1])
    filename = sys.argv[2]
    positions_x = [[0.4, 0.0, 0.7],[0.4, 0.35, 0.7], [0.5, 0.1, 0.5]]
    positions_y = [[0.5, 0.15, 0.6], [0.1, 0.7, 0.3], [0.7, 0.2, 0.9]]
    
    #Read all precision and recall and save to a hexbin plot.
    num_plots = 3
    for i in range(num_plots):
        [precision_weak, recall_weak, precision_promoter, recall_promoter, precision_enhancer, recall_enhancer] = read_all_precision_and_recall(cell[i], iters)
        plt.subplot(2, 2, i+1)
        make_hexbin_plots(precision_weak, recall_weak, precision_promoter, recall_promoter, precision_enhancer, recall_enhancer, cell[i], positions_x[i], positions_y[i], i)
    #Save the plots.
    plt.savefig(filename + ".png")
    
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
        directory = "/fs/project/PAS0272/Tara/DNase_SOM/" + cell + "/precision_recall_" + str(i)
        
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
def make_hexbin_plots(pw, rw, pp, rp, pe, re, cell, xpos, ypos, index):
    
    #Make plots for each regulatory category.
    sns.kdeplot(pw, rw, cmap=plt.cm.Greys)
    nonzeros = (np.nonzero(pe)[0]).tolist()
    zeros = (np.where(np.asarray(pe) == 0.0))
    try:
        sns.kdeplot(pe, re, cmap=plt.cm.Greys)
    except:
        sns.kdeplot(np.asarray(pe)[nonzeros], np.asarray(re)[nonzeros], cmap=plt.cm.Greys)
        #sns.regplot(x=np.array([0]), y=np.array([0]), scatter=True, fit_reg=False, marker='o',
        #    scatter_kws={"s": 10}, color = "black") 
    sns.kdeplot(pp, rp, cmap=plt.cm.Greys)
    plt.text(xpos[0], ypos[0], "Promoter")
    plt.text(xpos[1], ypos[1], "Enhancer")
    plt.text(xpos[2], ypos[2], "Weak")
    
    #Set x and y axis names.
    if index == 2:
        #Only show precision and recall for the bottom plot.
        plt.xlabel("Precision")
        plt.ylabel("Recall")
        #Add 'C' label.
        plt.text(-0.15, 1.1, "C", size=18)
    elif index == 0:
        #Remove x axis ticks for first plot.
        plt.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False) 
        #Add 'A' label.
        plt.text(-0.2, 1.1, "A", size=18)
    elif index == 1:
        #Remove y axis ticks for second plot.
        plt.tick_params(
            axis='y',
            which='both',
            left=False,
            right=False,
            labelleft=False) 
        #Make sure ticks match other plots.
        plt.xticks(np.arange(0, 1.5, 0.5))
        #Add 'B' label.
        plt.text(-0.1, 1.1, "B", size=18)
    plt.title(cell)

    
if __name__ == "__main__":
    main()