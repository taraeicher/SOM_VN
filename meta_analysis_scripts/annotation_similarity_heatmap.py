import numpy as np
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

"""
Output a heatmap of the counts of each type of match.
"""
def main():

    #Files needed for program
    all_shapes_path = sys.argv[1]
    match_log_path = sys.argv[2]
    ratio_out_path = sys.argv[3]
    heatmap_out_path = sys.argv[4]
    merged_path = sys.argv[5]
    exclude_enhancer = int(sys.argv[6]) > 0
    cell = sys.argv[7]

    #Compute the count of each annotation and annotation match.
    annotations = ["Promoter", "Enhancer", "Weak", "Unknown"]
    if exclude_enhancer:
        annotations = ["Promoter", "Weak", "Unknown"]
    count_each = get_total_counts(all_shapes_path, annotations, 0)
    count_merged = get_total_counts(merged_path, annotations, 1)
    match_counts = get_all_counts(match_log_path, annotations, count_each)
        
    #Print message to user.
    plot_heatmap(match_counts, annotations, heatmap_out_path, cell)
    
#Get the total count of clusters with each annotation.
def get_total_counts(file_path, names, anno_col):
    
    #Initialize counts and open file.
    counts = np.zeros(len(names))
    file = open(file_path, "r")
    #Count annotations.
    next_line = file.readline()
    while next_line:
        next_annotation = next_line.split("\t")[anno_col]
        try:
            i = names.index(next_annotation)
            counts[i] = counts[i] + 1
        except:
            pass
        next_line = file.readline()
    return counts
    
#Get the total count of matches with each annotation.
def get_all_counts(log_path, names, total_counts):
    
    #Initialize counts and open file.
    counts = np.zeros((len(names),len(names)))
    file = open(log_path, "r")
    next_line = file.readline()
    prev_name_1 = ""
    while next_line:
        
        #Get each component needed from the line.
        #The format is: Match between Weak and Weak(Brain_4_0) (Brain_13_0)
        split_line = next_line.split(" ")
        anno_1 = split_line[2]
        anno_2_and_name = next_line[4]
        anno_2 = split_line[4].split("(")[0]
        name_1 = split_line[4].split("(")[1].split(")")[0]
        name_2 = split_line[5][split_line[5].find("(")+1:split_line[5].find(")")]
        
        #Increment the positions.
        i = names.index(anno_1)
        j = names.index(anno_2)
        counts[i,j] = counts[i,j] + 1
        counts[j,i] = counts[j,i] + 1
        
        #Get the next line.
        next_line = file.readline()
    
    #Return counts.
    return counts   

#Plot a heatmap of total counts and a stacked bar plot of annotation type.
def plot_heatmap(matches, names, heat_path, cell):
    
    #Plot a heatmap of total counts.
    f = plt.figure(1, figsize=(8, 5))
    all_totals = np.transpose(np.tile(np.sum(matches, 0), len(matches)).reshape((len(matches), len(matches))))
    matches_scaled = matches / all_totals
    ax = plt.axes()
    heatmap = sns.heatmap(matches_scaled, cbar=True, cmap="binary", fmt="d", vmin = 0, vmax = 1, xticklabels = names, yticklabels = names, cbar_kws={'label': 'Percentage Cross-Correlated'})
    ax.set_title(cell)
    plt.xlabel("Regulatory Category of Cross-Correlated Shapes")
    plt.ylabel("Regulatory Category of Shape")
    plt.yticks(rotation = 0)
    fig = heatmap.get_figure()
    fig.savefig(heat_path + ".png")

if __name__ == "__main__":
    main()