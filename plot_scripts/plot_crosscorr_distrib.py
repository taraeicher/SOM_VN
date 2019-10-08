"""
Use the BED files containing each region's correlation to its closest shape.
Plot the correlation across all chromosomes, for the vanilla SOM, the permuted
WIG SOM, and the VNSSOM. The following arguments are required:
1. VNSSOM BED file directory containing annotations
2. Permuted SOM BED file directory containing annotations
3. Vanilla SOM BED file directory containing annotations
4. File where output plot should be saved.
"""

#Import sys for obtaining command line args.
import sys
import numpy as np

#Import matplotlib for plots.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    #Window size list and list of chromosomes are variable.
    #Allow user to choose number of iterations, initial learning
    #rate, and initial neighborhood size.
    #20, 50, 100, 250, 
    input_dir = sys.argv[1]
    input_dir_rand = sys.argv[2]
    input_dir_som = sys.argv[3]
    output = sys.argv[4]
    
    #Get all cross-correlations.
    chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    cc_vnssom = []
    cc_rand_vnssom = []
    cc_som = []
    for c in chroms:
        
        #Get all cross-corr for the given annotation in that file.
        in_file = open(input_dir + "/" + c, 'r')
        in_file_rand = open(input_dir_rand + "/" + c, 'r')
        in_file_som = open(input_dir_som + "/" + c, 'r')
        cc_vnssom.extend(get_crosscorrs(in_file))
        cc_rand_vnssom.extend(get_crosscorrs(in_file_rand))
        cc_som.extend(get_crosscorrs(in_file_som))
        in_file.close()
        in_file_rand.close()
        in_file_som.close()

    #Plot all cross-correlations.
    plot_densities(cc_vnssom, cc_rand_vnssom, cc_som, output)

"""
Extract all cross-correlations from a file.
"""
def get_crosscorrs(in_file):
    
    #Read in each line in the file and add to the cross-correlation list if it is the correct annotation.
    cross_corrs = []
    next_line = in_file.readline()
    while(next_line):
    
        #Get the chromosome and position info to print to the output file.
        #Get the input signal to match with the clusters.
        split_line = next_line.split("\t")
        cross_corrs.append(float(split_line[4]))           
        next_line = in_file.readline()
        
    return cross_corrs

"""
Plot the density distribution of cross-correlations for both the actual clusters and the randomized clusters.
"""
def plot_densities(vnssom, rand, som, output, cell):

        vnssom_plot = sns.distplot(vnssom, color = "black", label = "Promoters - Real WIG")
        rand_plot = sns.distplot(rand, color = "gray", label = "Enhancers - Real WIG")
        som_plot = sns.distplot(som, color = "gray", hist = False, kde_kws={'linestyle':'--'}, label = "Weak - Real WIG")
   
    plt.xlabel("Cross-Correlation")
    plt.ylabel("Density")
    fig = plotweak2.get_figure()
    fig.savefig(output)
    plt.close(fig)

if __name__ == "__main__":
    main()