import numpy as np
import pandas as pd
import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import glob, os

"""
Get all intensities from each WIG file and create a density plot of them.
"""
def main():

    #Get input data.
    a549_wig = sys.argv[1]
    h1_wig = sys.argv[2]
    brain_wig = sys.argv[3]
    out_file = sys.argv[4]
    
    #Get all intensity values.
    a549 = get_all_intensities(a549_wig)
    h1 = get_all_intensities(h1_wig)
    brain = get_all_intensities(brain_wig)
    
    #Plot densities.
    plot_densities(a549, h1, brain, out_file)
    
"""
Return a numeric list of all intensities for a WIG file.
"""
def get_all_intensities(wig_dir):
    
    #Loop through all chromosomes in the directory.
    full_list = []
    for wig_f in glob.glob(wig_dir + "*"):
    
        #Open the file.
        wig = open(wig_f, "r")
        
        #Discard the header lines.
        junk = wig.readline()
        junk = wig.readline()
        
        #Read all lines.
        lines = []
        line = wig.readline()
        while line:
            #Append the intensity value for each line.
            float_val = float(line.split("\t")[1])
            if float_val > 0:
                lines.append(math.log10(float_val))
            #Get next line.
            line = wig.readline()
        
        #Append the lines from this file to the full list.
        full_list.extend(lines)
        
    #Return the full list.
    return np.asarray(full_list)

"""
Plot the density distribution of cross-correlations for both the actual clusters and the randomized clusters.
"""
def plot_densities(a549, h1, brain, output):
    
    #Create the data frame.
    plt.figure()
    plt.title("Distribution of RPKM Signal Per Cell Type")
    rpkm_all = np.concatenate((a549, brain, h1))
    categories = np.concatenate((np.tile("A549", len(a549)), np.tile("Brain", len(brain)), np.tile("H1", len(h1))))
    df = pd.DataFrame({"Cell Type": categories, "Averaged RPKM Intensity Over 50 bp": rpkm_all})
    
    #Make and save violin plot containing all data.
    sns.violinplot(x="Cell Type", y="Averaged RPKM Intensity Over 50 bp", data=df, color="white")
    plt.savefig(output + "all.png") 
    
    #Make and save violin plot containing all data.
    plt.figure()
    plt.title("Distribution of RPKM Signal Per Cell Type")
    #df_50 = df[df["Averaged RPKM Intensity Over 50 bp"] <= 10]  
    #df_0 = df_50[df_50["Averaged RPKM Intensity Over 50 bp"] > 0]  
    sns.violinplot(x="Cell Type", y="Averaged RPKM Intensity Over 50 bp", data=df, color="white")
    plt.savefig(output + "wig_distrib.png") 
    
if __name__ == "__main__":
    main()