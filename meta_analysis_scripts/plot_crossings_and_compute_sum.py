"""
For each parameterization of the signal splitting process, plot
the distribution of cross counts. Create a final plot of the
sums of all crossings by parameterization.

The following parameters are required:
1. The input directory
2. The file name for the plot of sums
3. The output directory
"""

import numpy as np
import sys
import os
import glob
import pandas as pd
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
from mpl_toolkits.mplot3d import Axes3D

def main():

    # command-line arguments
    in_dir = sys.argv[1]
    plot_file_name = sys.argv[2]
    out_dir = sys.argv[3]
    
    # Create the output directory if it does not exist.
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # For each parameterization folder, read the crossings
    # files by chromosome. Plot and compute the sum.
    sums = dict()
    all_directories = glob.glob(in_dir + "/*_crossings*")
    for directory in tqdm(all_directories):
        crossings = []
        for file in glob.glob(directory + "/*"):
            crossings.extend(pd.read_csv(file, header = None).loc[:,0].values.tolist())
        sums[directory] = compute_sum(crossings)
        plot_crossings(crossings, directory, out_dir)
    plot_sums(sums, plot_file_name)
        
"""
Sum together crossing counts over all chromosomes.
"""
def compute_sum(crossings):
    crossings_array = np.array(crossings)
    return(np.sum(crossings_array))
  
"""
Create scatterplots of all sums by factor and by 
""" 
def plot_sums(sums, out_file):

    # Build a data frame of the sums by factor and margin.
    sums_df = np.zeros((len(sums.values()), 3))
    files = list(sums.keys())
    for i in range(len(files)):
        sums_df[i, 2] = sums[files[i]]
        local_name = split_file_name(files[i])
        sums_df[i, 0] = float(local_name.split("_")[0])
        sums_df[i, 1] = float(local_name.split("_")[1])
        
    # Make the scatterplot.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(sums_df[:,1], sums_df[:,0], sums_df[:,2], out_file)
    
    # Make the dimensions, etc, "pretty".
    ax.set_xlabel("Margin", labelpad = 15)
    ax.set_ylabel("Factor", labelpad = 15)
    ax.set_zlabel("Sum of Crossings", labelpad = 25)
    ax.tick_params(axis='z', which='major', pad=15)
    fig.set_size_inches(15,7, forward=True)
    
    # Save the plot.
    plt.savefig(out_file, bbox_inches = "tight")
    plt.close()
  
"""
Create a histogram of crossings across all chromosomes in a directory.
""" 
def plot_crossings(crossings, dir, out_dir):

    local_name = split_file_name(dir)
    plt.hist(crossings, color = "black", bins = 100)
    plt.xlim(0, 50)
    plt.ylim(0, 150000)
    plt.savefig(out_dir + "/" + local_name + ".png", bbox_inches = "tight")
    plt.close()
    
"""
Split the file name and return the local piece (factor + margin)
"""
def split_file_name(dir):
    split_by_slash = dir.split("/")
    split_by_underscore = split_by_slash[len(split_by_slash) - 1].split("_")
    return(split_by_underscore[0] + "_" + split_by_underscore[1])
    
if __name__ == "__main__":
    main()