#Import sys for obtaining command line args.
import sys
import numpy as np

#Import matplotlib for plots.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
matplotlib.rcParams.update({'font.size': 24})

"""
For each datum in each of the window size files, match it to its
closest cluster in the cluster files. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

    #Window size list and list of chromosomes are variable.
    #Allow user to choose number of iterations, initial learning
    #rate, and initial neighborhood size.
    #20, 50, 100, 250,
    input_dir = sys.argv[1]
    input_dir_rand = sys.argv[2]
    output = sys.argv[3]
    cell = sys.argv[4]

    #Get all cross-correlations.
    chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    cc_prom = []
    cc_enh = []
    cc_weak = []
    #rand_cc_prom = []
    #rand_cc_enh = []
    rand_cc_weak = []
    for c in chroms:

        #Get all cross-corr for the given annotation in that file.
        in_file = open(input_dir + c, 'r')
        in_file_rand = open(input_dir_rand + c, 'r')
        cc_prom.extend(get_crosscorr_anno(in_file, "Promoter"))
        cc_enh.extend(get_crosscorr_anno(in_file, "Enhancer"))
        cc_weak.extend(get_crosscorr_anno(in_file, "Weak"))
        #rand_cc_prom.extend(get_crosscorr_anno(in_file_rand,"Promoter"))
        #rand_cc_enh.extend(get_crosscorr_anno(in_file_rand,"Enhancer"))
        rand_cc_weak.extend(get_crosscorr_anno(in_file_rand,"Weak"))
        in_file.close()
        in_file_rand.close()
   
    #Plot all cross-correlations.
    plot_densities(cc_prom, cc_enh, cc_weak, rand_cc_weak, output, cell)

"""
Match each input with the given window size to the nearest cluster for that window size.
Print out the region with its corresponding cluster to a BED file.
"""
def get_crosscorr_anno(in_file, annotation):

    #Read in each line in the file and add to the cross-correlation list if it is the correct annotation.
    cross_corrs = []
    next_line = in_file.readline()
    while(next_line):

        #Get the chromosome and position info to print to the output file.
        #Get the input signal to match with the clusters.
        split_line = next_line.split("\t")
        anno = split_line[3]
        cross_corr = float(split_line[4])

        #Match the data to the nearest cluster and obtain the match and the ambiguity metric.
        if anno == annotation:
            cross_corrs.append(cross_corr)

        #Read the next line in the file.
        next_line = in_file.readline()

    #Go bac to the beginning of the file and return the cross-corr list.
    in_file.seek(0)
    return cross_corrs

"""
Plot the density distribution of cross-correlations for both the actual clusters and the randomized clusters.
"""
def plot_densities(promoters, enhancers, weaks, rand_weaks, output, cell):

    #Plot promoters, enhancers, weak, and permuted weak.
    plotprom1 = sns.distplot(promoters, color = "green")
    plotenh1 = sns.distplot(enhancers, color = "blue")
    plotweak1 = sns.distplot(weaks, color = "red")
    plotweak2 = sns.distplot(rand_weaks, color = "black", kde_kws={'linestyle':'--'})
    plt.title(cell)
    plt.xlabel("Cross-Correlation")
    plt.ylabel("Density")
    fig = plotweak2.get_figure()
    fig.savefig(output + ".png", bbox_inches = "tight", dpi = 300)
    plt.close(fig)

    legend_elements = [Line2D([0], [0],color="blue", label='Enhancer (Original Signal)'),
                        Line2D([0], [0],color="green", label='Promoter (Original Signal)'),
                        Line2D([0], [0],color="red", label='Weak (Original Signal)'),
                        Line2D([0], [0], color = "black", linestyle="--", label='Weak (Permuted Signal)')]
    plt.axis("off")
    plt.legend(handles=legend_elements, loc="lower right")

    #Save the plot.
    plt.savefig(output + "_legend.png", bbox_inches = "tight", dpi = 300)

if __name__ == "__main__":
    main()
