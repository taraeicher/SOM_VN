import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size': 20})


"""
For each of the annotations, find its information gain for each cluster.
For those cluster-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the bed file and cluster data.
    bed_dir = sys.argv[1]
    wig_dir = sys.argv[2]
    sig_dir = sys.argv[3]
    unknown_dir = sys.argv[4]
    output_dir = sys.argv[5]
    training = sys.argv[6]
    cell = sys.argv[7]
    chrom_names = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', 'X', 'Y']
    predictions = ["Promoter", "Enhancer", "Weak", "Other", "Unknown"]
    
    #Save unknowns.
    save_unknown_chart(unknown_dir, output_dir, cell, chrom_names, training)
    
#Save a line with the ChromHMM distribution for each unknown annotation.  
def save_unknown_chart(file_dir, output, cell, chrom_names, training):
        
    #Get all percentages for each chromosome.
    unknown_percentage = []
    promoter_percentage = []
    enhancer_percentage = []
    weak_percentage = []
    polycomb_percentage = []
    other_percentage = []

    for c in chrom_names:
        f = open(file_dir + "/unknown_dist_" + c, "r")
        unknown_percentage.append(f.readline())
        promoter_percentage.append(f.readline())
        enhancer_percentage.append(f.readline())
        weak_percentage.append(f.readline())
        other_percentage.append(f.readline())
        f.close()
    unknowns = np.asarray([float(i) for i in unknown_percentage])
    promoters = np.asarray([float(i) for i in promoter_percentage])
    enhancers = np.asarray([float(i) for i in enhancer_percentage])
    weaks = np.asarray([float(i) for i in weak_percentage])
    other = np.asarray([float(i) for i in other_percentage])

    #Plot each line.
    r = range(0, len(chrom_names))
    plt.figure(figsize=(12,6)) 
    # Create promoter line
    plt.plot(r, promoters, color = "black", label = "Promoter")
    # Create enhancer line
    plt.plot(r,  enhancers, color = "gray", label = "Enhancer")
    # Create weak line
    plt.plot(r, weaks, color = "silver", label = "Weak")
    # Create other line
    plt.plot(r, other, color = "gray", linestyle = '--', label = "Other")
    
    # Add title, axis, and legend. Save plot.
    plt.xticks(r, chrom_names)
    plt.xlabel("Chromosome")
    plt.ylabel("Percentage")
    plt.title(training + " to " + cell)
    plt.savefig(output + "/unknown_" + training + ".png")
    plt.close()
    
    legend_elements = [Line2D([0], [0], marker=None, label='Enhancer', color = 'gray'),
                        Line2D([0], [0], marker=None, label='Promoter', color = 'black'),
                        Line2D([0], [0], marker=None, label='Weak', color = 'silver'),
                        Line2D([0], [0], marker=None, label='Other', color = 'gray', linestyle = '--'),
                       ]
    plt.legend(handles=legend_elements, loc="lower right")
    plt.savefig(output + "/legend.png")
    plt.close()
    
if __name__ == "__main__":
    main()