import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
from matplotlib import rc
import glob
import math

"""
For each of the annotations, find its information gain for each cluster.
For those cluster-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
Y_MAX = 50
WINDOW_SIZE = 80
def main():
    
    #Get list of clusters.
    clusters = []
    cluster_file = open(sys.argv[1], 'r')
    cell = sys.argv[9]
    next_clust = cluster_file.readline()
    while next_clust:
        clusters.append(next_clust)
        next_clust = cluster_file.readline()
    
    #Get distribution of cluster assignment and TSS sites per cluster.
    #Then, plot the clusters.
    annotation_files_brain = sorted(glob.glob(sys.argv[3] + "*" + "clust.bed"))
    annotation_files_a549 = sorted(glob.glob(sys.argv[4] + "*" + "clust.bed"))
    annotation_files_h1 = sorted(glob.glob(sys.argv[5] + "*" + "clust.bed"))
    tss_files_brain = sorted(glob.glob(sys.argv[6] + "*" + ".bed"))
    tss_files_a549 = sorted(glob.glob(sys.argv[7] + "*" + ".bed"))
    tss_files_h1 = sorted(glob.glob(sys.argv[8] + "*" + ".bed"))
    anno_counts_brain = get_annotation_distribution(clusters, annotation_files_brain)
    anno_counts_a549 = get_annotation_distribution(clusters, annotation_files_a549)
    anno_counts_h1 = get_annotation_distribution(clusters, annotation_files_h1)
    tss_avg_brain = get_expected_tss(clusters, tss_files_brain, anno_counts_brain)
    tss_avg_a549 = get_expected_tss(clusters, tss_files_a549, anno_counts_a549)
    tss_avg_h1 = get_expected_tss(clusters, tss_files_h1, anno_counts_h1)
    save_line_charts(clusters, sys.argv[2], anno_counts_brain, anno_counts_a549, anno_counts_h1, tss_avg_brain, tss_avg_a549, tss_avg_h1, cell)

#Create a dictionary with count of cluster annotations for each cluster.
def get_annotation_distribution(clusters, annotated_bed_list):
    cluster_counts = dict()

    #For each line in the file, track its annotation.
    for annotated_bed in annotated_bed_list:
        bed = np.genfromtxt(annotated_bed, delimiter = "\t", dtype = str)
        for line in bed:
            anno = line[3]
            if anno in cluster_counts:
                cluster_counts[anno] += 1
            else:
                cluster_counts[anno] = 1
   
    #Transform the counts into percentages.
    #cluster_perc = {k: v / total for k, v in cluster_counts.iteritems()}
    return cluster_counts
    
#Create a dictionary with percent of clusters containing TSS.
def get_expected_tss(clusters, tss_bed_list, anno_counts):
    cluster_counts = dict()
    
    #For each line in the file, track its annotation.
    for tss_bed in tss_bed_list:
        bed = np.genfromtxt(tss_bed, delimiter = "\t", dtype = str)
        for line in bed:
            anno = line[3]
            tss = int(line[5])
            if anno in cluster_counts and tss >= 1:
                cluster_counts[anno] += tss
            elif not(anno in cluster_counts):
                cluster_counts[anno] = tss
    
    #Transform the counts into percentages.
    cluster_perc = dict()
    for anno in cluster_counts:
        cluster_perc[anno] = cluster_counts[anno] / anno_counts[anno]
    return cluster_perc

#Save the line chart of each important cluster.  
def save_line_charts(clusters, path, anno_brain, anno_a549, anno_h1, tss_brain, tss_a549, tss_h1, cell):
   
    #Set up colors.
    colors = dict()
    colors["Promoter"] = 'black'
    colors["Enhancer"] = 'gray'
    colors["Weak"] = 'silver'
    
    #Set up initial plots.
    r = range(0, WINDOW_SIZE)
    for i in range(1, 8):
        plt.figure(i, figsize=(12,6)) 
        axes = plt.gca()
        axes.set_xlim([0,WINDOW_SIZE])
        axes.set_ylim([0,Y_MAX])
        plt.xlabel("Position in Shape (" + str(BIN_SIZE) + " bp)")
        plt.ylabel("RPKM Intensity")
        plt.title("Promoters - " + str(cell))
    
    for i in range(8, 11):
        plt.figure(i, figsize=(12,6)) 
        axes = plt.gca()
        axes.set_xlim([0,WINDOW_SIZE])
        axes.set_ylim([0,Y_MAX])
        plt.xlabel("Position in Shape (" + str(BIN_SIZE) + " bp)")
        plt.ylabel("RPKM Intensity")
        plt.title("Enhancers - " + str(cell))
    
    plt.figure(11, figsize=(12,6)) 
    axes = plt.gca()
    axes.set_xlim([0,WINDOW_SIZE])
    axes.set_ylim([0,Y_MAX])
    plt.xlabel("Position in Shape (" + str(BIN_SIZE) + " bp)")
    plt.ylabel("RPKM Intensity")
    plt.title("Weak - " + str(cell))
        
    #Plot all clusters.
    i = 0
    j = 0
    for c in range(0, len(clusters)):
        pieces = clusters[c].split("\t")
        name = pieces[0]
        annotation = pieces[1]
        signal = np.array([float(i) for i in pieces[2].split(",")])
        
        #Plot everything that isn't unknown, which the corresponding color
        #to match the other figures.
        
        if annotation == "Promoter":
            
            #Plot the figure.
            i += 1
            plt.figure(math.ceil(i / 10), figsize=(12,6))
            text_pos_x = np.argmax(signal)
            text_pos_y = min(np.max(signal), Y_MAX - 2)
            plt.plot(r, signal, color = colors[annotation], linewidth = 4)
            
            #Change the plot coordinates for specific clusters.
            if name == "A549_17_7":
                text_pos_y = 43.0
            elif name == "A549_17_3":
                text_pos_x = 20.0
            elif name == "A549_17_5":
                text_pos_x = 20.0
                text_pos_y = 40.0
            elif name == "A549_17_8":
                text_pos_y = 30.0
            elif name == "A549_19_18":
                text_pos_x = 25.0
            elif name == "A549_19_10":
                text_pos_x = 30.0
                text_pos_y = 30.0
            elif name == "A549_21_0":
                text_pos_y = 40.0
            elif name == "A549_17_38":
                text_pos_y = 5.0
            elif name == "A549_17_20":
                text_pos_x = 23.0
                text_pos_y = 23.0
            elif name == "A549_17_22":
                text_pos_y = 20.0
            elif name == "Brain_17_10":
                text_pos_y = 25.0
            elif name == "A549_17_29":
                text_pos_y = 28.0
            elif name == "A549_19_13":
                text_pos_x = 53.0
            elif name == "A549_17_1":
                text_pos_y = 8.0
            elif name == "A549_18_15":
                text_pos_y = 35.0
            elif name == "A549_20_23":
                text_pos_y = 22.0
            elif name == "A549_14_28":
                text_pos_y = 5.0
            elif name == "A549_20_29":
                text_pos_x = 45.0
            elif name == "A549_9_2":
                text_pos_y = 40.0
            elif name == "A549_6_5":
                text_pos_y = 10.0
            elif name == "A549_12_3":
                text_pos_y = 40.0
            elif name == "Brain_21_10":
                text_pos_y = 38.0
            elif name == "Brain_21_5":
                text_pos_y = 43.0
            elif name == "Brain_20_2":
                text_pos_y = 30.0
                text_pos_x = 47.0
            elif name == "Brain_20_0":
                text_pos_y = 42.5
                text_pos_x = 20.0
            elif name == "Brain_21_2":
                text_pos_y = 38.0
                text_pos_x = 20.0
            elif name == "Brain_12_24":
                text_pos_y = 35.0
            elif name == "Brain_12_1":
                text_pos_y = 42.0
            elif name == "Brain_12_4":
                text_pos_x = 22.0
            elif name == "Brain_3_2":
                text_pos_y = 35.0
            elif name == "Brain_21_7":
                text_pos_y = 42.0
                text_pos_x = 25.0
            elif name == "H1_22_3":
                text_pos_y = 24.0
            elif name == "H1_3_2":
                text_pos_y = 40.0
                text_pos_x = 0.0
            elif name == "H1_21_3":
                text_pos_x = 60.0
            elif name == "H1_19_3":
                text_pos_x = 60.0
            elif name == "H1_16_8":
                text_pos_x = 70.0
            elif name == "H1_19_6":
                text_pos_x = 20.0
                text_pos_y = 7.0
            elif name == "H1_16_2":
                text_pos_y = 42.0
                
            #Add TSS and count information.
            plt.text(text_pos_x, text_pos_y, name, weight = "bold")
            try:
                plt.text(text_pos_x, text_pos_y - 1.25, "Brain Count: " + str(anno_brain[name]) + ", TSS: " + "{:.2f}".format(tss_brain[name]))
                plt.text(text_pos_x, text_pos_y - 2.5, "A549 Count: " + str(anno_a549[name]) + ", TSS: " + "{:.2f}".format(tss_a549[name]))
                plt.text(text_pos_x, text_pos_y - 3.75, "H1 Count: " + str(anno_h1[name]) + ", TSS: " + "{:.2f}".format(tss_h1[name]))
            except:
                pass
                
        elif annotation == "Enhancer":
            j += 1
            #Plot the figure.
            plt.figure(math.ceil(j / 10) + 7, figsize=(12,6))
            text_pos_x = np.argmax(signal)
            text_pos_y = min(np.max(signal), Y_MAX - 2)
            plt.plot(r, signal, color = colors[annotation], linewidth = 4)
            
            #Change the plot coordinates for specific clusters.
            if name == "A549_20_34":
                text_pos_y = 12.0
            elif name == "Brain_20_18":
                text_pos_y = 50.0
            elif name == "Brain_19_18":
                text_pos_x = 60.0
                text_pos_y = 10.0
            elif name == "Brain_20_15":
                text_pos_x = 60.0
            elif name == "Brain_20_8":
                text_pos_x = 23.0
            elif name == "Brain_20_2":
                text_pos_x = 42.0
                text_pos_y = 25.0
            elif name == "Brain_19_18":
                text_pos_y = 17.0
                text_pos_x = 43.0
            elif name == "Brain_19_7":
                text_pos_x = 42.0
                text_pos_y = 10.0
            elif name == "Brain_19_27":
                text_pos_x = 60.0
                text_pos_y = 15.0
            elif name == "Brain_19_4":
                text_pos_y = 5.0
                text_pos_x = 1.0
            elif name == "Brain_19_14":
                text_pos_y = 15.0
                text_pos_x = 42.0
            elif name == "Brain_17_25":
                text_pos_x = 20.0
            elif name == "Brain_17_5":
                text_pos_x = 45.0
                text_pos_y = 13.0
            elif name == "Brain_17_17":
                text_pos_y = 30.0
            elif name == "Brain_17_20":
                text_pos_y = 23.0
            elif name == "Brain_11_9":
                text_pos_y = 12.0
                text_pos_x = 45.0
            elif name == "Brain_12_1":
                text_pos_y = 40.0
            elif name == "Brain_12_4":
                text_pos_x = 22.0
            elif name == "Brain_16_0":
                text_pos_y = 18.0
            elif name == "Brain_16_12":
                text_pos_x = 25.0
                text_pos_y = 12.0
                
            #Add the TSS and count information.
            plt.text(text_pos_x, text_pos_y, name, weight = "bold")
            try:
                plt.text(text_pos_x, text_pos_y - 1.25, "Brain Count: " + str(anno_brain[name]) + ", TSS: " + "{:.2f}".format(tss_brain[name]))
                plt.text(text_pos_x, text_pos_y - 2.5, "A549 Count: " + str(anno_a549[name]) + ", TSS: " + "{:.2f}".format(tss_a549[name]))
                plt.text(text_pos_x, text_pos_y - 3.75, "H1 Count: " + str(anno_h1[name]) + ", TSS: " + "{:.2f}".format(tss_h1[name]))
            except:
                pass            
            
        elif annotation == "Weak":
            
            #Plot the figure.
            plt.figure(11, figsize=(12,6))
            text_pos_x = np.argmax(signal)
            text_pos_y = min(np.max(signal), Y_MAX - 2)
            plt.plot(r, signal, color = colors[annotation], linewidth = 4)
            
            #Change the plot coordinates for specific clusters.
            if name == "A549_19_18":
                text_pos_x = 60.0
                text_pos_y = 5.0
            elif name == "Brain_13_16":
                text_pos_x = 20.0
                test_pos_y = 10.0
            elif name == "Brain_13_5":
                text_pos_x = 1.0
                text_pos_y = 10.0
            elif name == "Brain_21_38":
                text_pos_x = 25.0
                text_pos_y = 10.0
            elif name == "Brain_21_20":
                text_pos_x = 60.0
                text_pos_y = 10.0
            elif name == "Brain_4_18":
                text_pos_y = 15.0
            elif name == "Brain_7_18":
                text_pos_y = 5.0
            elif name == "Brain_4_10":
                text_pos_x = 60.0
                text_pos_y = 15.0
            elif name == "Brain_3_17":
                text_pos_x = 60.0
            elif name == "Brain_21_16":
                text_pos_y = 10.0
            elif name == "Brain_21_18":
                text_pos_y = 5.0
            elif name == "Brain_4_19":
                text_pos_y = 17.0
            elif name == "Brain_19_28":
                text_pos_x = 20.0
                text_pos_y = 12.0
            elif name == "H1_20_15":
                text_pos_x = 62.0
                text_pos_y = 5.0
            elif name == "H1_10_3":
                text_pos_y = 5.0
            elif name == "H1_21_6":
                text_pos_x = 20.0
                text_pos_y = 13.0
            elif name == "H1_18_6":
                text_pos_x = 0.0
                text_pos_y = 13.0
            elif name == "H1_9_14":
                text_pos_x = 0.0
                text_pos_y = 5.0
            elif name == "H1_7_13":
                text_pos_x = 20.0
                text_pos_y = 5.0
            
            #Add TSS and count information.
            plt.text(text_pos_x, text_pos_y, name, weight = "bold")
            try:
                plt.text(text_pos_x, text_pos_y - 1.25, "Brain Count: " + str(anno_brain[name]) + ", TSS: " + "{:.2f}".format(tss_brain[name]))
                plt.text(text_pos_x, text_pos_y - 2.5, "A549 Count: " + str(anno_a549[name]) + ", TSS: " + "{:.2f}".format(tss_a549[name]))
                plt.text(text_pos_x, text_pos_y - 3.75, "H1 Count: " + str(anno_h1[name]) + ", TSS: " + "{:.2f}".format(tss_h1[name]))
            except:
                pass

    #Save all plots. 
    for f in range(1, 8):
        plt.figure(f, figsize=(12,6))
        plt.savefig(path + str(cell) + "_Promoters" + str(f) + ".png")
        plt.close()
    for f in range(8, 11):
        plt.figure(f, figsize=(12,6))
        plt.savefig(path + str(cell) + "_Enhancers" + str(f) + ".png")
        plt.close()
    plt.figure(11, figsize=(12,6))
    plt.savefig(path + str(cell) + "_Weak.png")
    plt.close()
    
if __name__ == "__main__":
    main()