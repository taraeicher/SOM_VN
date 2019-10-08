import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import os
plt.rcParams.update({'font.size': 14})
sys.path.append(os.path.abspath("../common_scripts"))
import wig_and_signal_utils as wsu

"""
Read in all annotations.
"""
def main():

    #Read in all data.
    all_regions = dict()
    file_in = open(sys.argv[1], "r")
    next_line = file_in.readline()
    while next_line:
        split_line = next_line.split("\t")
        floats = np.array([float(s) for s in split_line[2].split(",")])
        if all_regions.get(split_line[0]) == None:
            all_regions[split_line[0]] = list()
        elif not np.isnan(floats).any(): 
            all_regions[split_line[0]].append(floats)
        next_line = file_in.readline()

    #Convert to numpy.
    for key in all_regions:
        all_regions[key] = np.stack(all_regions[key])

    #Plot everything.
    plot_bars(sys.argv[2], all_regions)
 
    
"""
Plot annotations in a stacked bar plot.
""" 
def plot_bars(out, percentages):

    #Plot each type of region separately.
    for key in percentages:
        
        #Set up parameters for plot.
        perc = percentages[key]
        bar_count = percentages[key].shape[0]
        fig, ax = plt.subplots(figsize=(5, 15))
        wid = 0.75
        color_polycomb = "white"
        color_promoter = "black"
        color_enhancer = "gray"
        color_weak = "white"
        
        #Plot shape-based predictions.
        tot_wid = 4000
        xs = np.arange(perc.shape[0])
        p1_s = plt.barh(xs, perc[:,2], color = color_polycomb, edgecolor = "black", hatch="+", height = wid)
        p3_s = plt.barh(xs, perc[:,0], color = color_promoter, edgecolor = "black", left = [i for i in perc[:,2]],  height = wid)
        p4_s = plt.barh(xs, perc[:,1], color = color_enhancer, edgecolor = "black", left = [i+j for i,j in zip(perc[:,2], perc[:,0])], height = wid)
        p5_s = plt.barh(xs, perc[:,3], color = color_weak, edgecolor = "black", left = [i+j+k for i,j,k in zip(perc[:,2], perc[:,0], perc[:,1])], height = wid)
        plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
        plt.box(on=None)
        
        # Save plot.
        plt.savefig(out + "_" + key + ".png")
        plt.close()
        
        #Save the legend.
        legend_elements = [Patch(facecolor=color_polycomb, edgecolor = "black", hatch = "+", label='Polycomb'),
                            Patch(facecolor=color_promoter, edgecolor = "black", label='Promoter'),
                            Patch(facecolor=color_enhancer, edgecolor = "black", label='Enhancer'),
                            Patch(facecolor=color_weak, edgecolor = "black", label='Weak')
                           ]
        plt.legend(handles=legend_elements, loc="lower right", title = "Percentage of Bins")
        plt.savefig(out + "_legend" + ".png")

    
if __name__ == "__main__":
    main()