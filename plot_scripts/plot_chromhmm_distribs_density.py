"""
Plot the percentage distributions of promoter, enhancer,
mnemonic, and weak in each region. Arguments:
1. ChromHMM distributions for GM12878
2. ChromHMM distributions for A549
3. ChromHMM distributions for Brain
4. ChromHMM distributions for H1
2. File where plots should be saved
"""

#Import sys for obtaining command line args.
import sys
import numpy as np
import pandas as pd

#Import matplotlib for plots.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    #Input and output files
    GM12878 = pd.read_csv(open(sys.argv[1], "r"))
    A549 = pd.read_csv(open(sys.argv[2], "r"))
    Brain = pd.read_csv(open(sys.argv[3], "r"))
    H1 = pd.read_csv(open(sys.argv[4], "r"))
    output = sys.argv[5]
    nonzero = sys.argv[6] == "True"

    #Plot all cross-correlations.
    plt.subplot(2,2,1)
    promoter_fig = plot_densities(GM12878, A549, Brain, H1, "Promoter", False, False, True, False, nonzero)
    plt.subplot(2,2,2)
    enhancer_plot = plot_densities(GM12878, A549, Brain, H1, "Enhancer", False, False, True, True, nonzero)
    plt.subplot(2,2,3)
    plot_densities(GM12878, A549, Brain, H1, "Repressor", True, True, True, False, nonzero)
    plt.subplot(2,2,4)
    plot_densities(GM12878, A549, Brain, H1, "Weak", False, True, True, False, nonzero)
    plt.savefig(output, bbox_inches = "tight", format = "pdf")


"""
Plot the percentage distribution for
each type of RE.
"""
def plot_densities(GM12878, A549, Brain, H1, RE, include_y, include_xtick, include_ytick, include_legend, nonzero):

    bin_locs=np.arange(0.0, 1.1, 0.1)
    threshold = -1
    if nonzero:
        threshold = 0
    plt.hist(GM12878.loc[:,RE].iloc[np.where(GM12878.loc[:,RE] > threshold)[0]], color = "red", label = "GM12878", alpha = 0.5, bins = bin_locs)
    plt.hist(A549.loc[:,RE].iloc[np.where(A549.loc[:,RE] > threshold)[0]], color = "orange", label = "A549", alpha = 0.5, bins = bin_locs)
    plt.hist(Brain.loc[:,RE].iloc[np.where(Brain.loc[:,RE] > threshold)[0]], color = "blue", label = "Brain", alpha = 0.5, bins = bin_locs)
    plt.hist(H1.loc[:,RE].iloc[np.where(H1.loc[:,RE] > threshold)[0]], color = "gray", label = "H1", alpha = 0.5, bins = bin_locs)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.xlabel("Percent " + RE)
    #plt.ylim(0, 2500000)
    if include_legend:
        plt.legend(bbox_to_anchor=(1.6, 1.05))
    if not include_ytick:
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    if include_y:
        plt.ylabel("Count of Regions")
    if not include_xtick:
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

if __name__ == "__main__":
    main()