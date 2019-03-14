import numpy as np
import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({'font.size': 18})
import seaborn as sns
import math

"""
Get all intensities from each WIG file and create a density plot of them.
"""
def main():

    #Get input data.
    a549_bed = open(sys.argv[1], "r")
    h1_bed = open(sys.argv[2], "r")
    brain_bed = open(sys.argv[3], "r")
    out_file = sys.argv[4]
    
    #Get all intensity values.
    a549_promoter, a549_enhancer, a549_weak = get_all_sizes(a549_bed)
    h1_promoter, h1_enhancer, h1_weak = get_all_sizes(h1_bed)
    brain_promoter, brain_enhancer, brain_weak = get_all_sizes(brain_bed)
    
    #Plot densities.
    plot_densities(a549_promoter, a549_enhancer, a549_weak, h1_promoter, h1_enhancer, h1_weak, brain_promoter, brain_enhancer, brain_weak, out_file)
    
"""
Return a numeric list of all sizes in the ChromHMM BED file.
"""
def get_all_sizes(bed):
    
    #Annotation-specific information to use in counting
    promoters = []
    enhancers = []
    weak = []
    promoter_annos = ["1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"]
    enhancer_annos = ["6_EnhG", "7_Enh", "12_EnhBiv"]
    weak_annos = ["9_Het", "15_Quies"]
    anno_col = 3
    
    #Loop through each line in the file and add its length to the appropriate list.
    line = bed.readline().split("\t")
    while len(line) > 1:
    
        #Compute the size of the region.
        size = int(line[2]) - int(line[1])
        anno = line[3].rstrip()
        
        #Add size to appropriate list.
        if anno in promoter_annos:
            promoters.append(math.log10(size))
        elif anno in enhancer_annos:
            enhancers.append(math.log10(size))
        elif anno in weak_annos:
            weak.append(math.log10(size))
            
        #Get next line.
        line = bed.readline().split("\t")
        
    #Return the full list.
    return [promoters, enhancers, weak]

"""
Plot the distribution of ChromHMM annotation sizes for all cell lines and annotation types.
"""
def plot_densities(a549_promoter, a549_enhancer, a549_weak, h1_promoter, h1_enhancer, h1_weak, brain_promoter, brain_enhancer, brain_weak, output):
    
    #Create the data frame.
    plt.figure(figsize=(8,8))
    sizes_all = np.concatenate((a549_promoter, a549_enhancer, a549_weak, brain_promoter, brain_enhancer, brain_weak, h1_promoter, h1_enhancer, h1_weak))
    categories = np.concatenate((np.tile("A549 (P)", len(a549_promoter)), np.tile("A549 (E)", len(a549_enhancer)), np.tile("A549 (W)", len(a549_weak)), np.tile("Brain (P)", len(brain_promoter)), np.tile("Brain (E)", len(brain_enhancer)), np.tile("Brain (W)", len(brain_weak)), np.tile("H1 (P)", len(h1_promoter)), np.tile("H1 (E)", len(h1_enhancer)), np.tile("H1 (W)", len(h1_weak))))
    types = np.concatenate((np.tile("Promoter", len(a549_promoter)), np.tile("Enhancer", len(a549_enhancer)), np.tile("Weak", len(a549_weak)), np.tile("Promoter", len(brain_promoter)), np.tile("Enhancer", len(brain_enhancer)), np.tile("Weak", len(brain_weak)), np.tile("Promoter", len(h1_promoter)), np.tile("Enhancer", len(h1_enhancer)), np.tile("Weak", len(h1_weak))))
    df = pd.DataFrame({"Cell Type": categories, "Log of Width in bp (Base 10)": sizes_all, "Regulatory Category": types})
    colors = ["black", "gray", "white", "black", "gray", "white", "black", "gray", "white"]
    
    #Make and save violin plot containing all data.
    #df_50 = df[df["Width in bp"] <= 250000]  
    sns.violinplot(x="Cell Type", y="Log of Width in bp (Base 10)", data = df, palette = colors, inner = None)
    plt.xticks((1,4,7), ('A549', 'Brain', 'H1'))
    #plt.xticks(rotation = 45)
    plt.savefig(output + ".png")
    plt.close()
    
    #Add the legend.
    legend_elements = [Patch(facecolor="gray", edgecolor="black", label='Enhancer'),
                        Patch(facecolor="black", edgecolor="black", label='Promoter'),
                        Patch(facecolor="white", edgecolor="black", label='Weak'),
                       ]
    #fig, ax = plt.subplots()
    plt.legend(handles=legend_elements, loc="lower right")
    
    #Save the plot.
    plt.savefig(output + "_legend.png")
    plt.close()
    
if __name__ == "__main__":
    main()