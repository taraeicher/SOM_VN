import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #For each chromosome, get all distributions.
    total_comp = 9
    percentages = np.zeros((9, 3, 3))
    cells = ["A549", "Brain", "H1"]
    idx = 0
    for src in cells:
        for dest in cells:
            if src != "H1":
                percentages[idx] = np.genfromtxt(sys.argv[1] + dest + "/" + "mispredictions" + src + ".csv", delimiter = ",")
            else:
                percentages[idx, (0,2)] = np.genfromtxt(sys.argv[1] + dest + "/" + "mispredictions" + src + ".csv", delimiter = ",")
            idx = idx + 1
    plot_out = sys.argv[2]
    labels = ["Promoter", "Enhancer", "Weak"]

    #Plot distributions in a bar plot for all cell combos.
    plot_bars(percentages, labels, plot_out, total_comp)
    
"""
Plot annotations in a stacked bar plot.
""" 
def plot_bars(data, labels, output, count):

    #Set up parameters for plot.
    bar_count = 4
    pos = list(range(0, count * bar_count, bar_count)) 
    fig, ax = plt.subplots(figsize=(10,6))
    wid = 0.75
    w = 1
    ax.set_xlim(-1, len(pos) * bar_count)
    ax.set_ylim(0, 1.08)
    color_promoter = "black"
    color_enhancer = "gray"
    color_weak = "white"
    names = ["A-A", "A-B", "A-H", "B-A", "B-B", "B-H", "H-A", "H-B", "H-H"]

    #Plot promoter ground truth.
    p1_s = plt.bar(pos, data[:,0,0], color = color_promoter, edgecolor = color_promoter, width = wid)
    p2_s = plt.bar(pos, data[:,0,1], color = color_enhancer, edgecolor = color_enhancer, bottom = data[:,0,0], width = wid)
    p3_s = plt.bar(pos, data[:,0,2], color = color_weak, edgecolor = "black", bottom = [i+j for i,j in zip(data[:,0,0], data[:,0,1])],  width = wid)
    for p in pos:
        ax.text(p, 1.0, "P", ha='center', va='bottom')
    
    #Plot enhancer ground truth.
    pos2 = [p + w for p in pos]
    p1_t = plt.bar(pos2, data[:,1,0], color = color_promoter, edgecolor = color_promoter, width = wid)
    p2_t = plt.bar(pos2, data[:,1,1], color = color_enhancer, edgecolor = color_enhancer, bottom = data[:,1,0], width = wid)
    p3_t = plt.bar(pos2, data[:,1,2], color = color_weak, edgecolor = "black", bottom = [i+j for i,j in zip(data[:,1,0], data[:,1,1])],  width = wid)
    for p in pos2:
        ax.text(p, 1.0, "E", ha='center', va='bottom')
        
    #Plot shape-based predictions AND position-based predictions.
    pos3 = [p + w * 2 for p in pos]
    p1_t = plt.bar(pos3, data[:,2,0], color = color_promoter, edgecolor = color_promoter, width = wid)
    p2_t = plt.bar(pos3, data[:,2,1], color = color_enhancer, edgecolor = color_enhancer, bottom = data[:,2,0], width = wid)
    p3_t = plt.bar(pos3, data[:,2,2], color = color_weak, edgecolor = "black", bottom = [i+j for i,j in zip(data[:,2,0], data[:,2,1])],  width = wid)
    for p in pos3:
        ax.text(p, 1.0, "W", ha='center', va='bottom')

    plt.xlabel("Training and Testing Pair")
    plt.ylabel("Percentage")
    #plt.legend(loc="center right")  

    # Add title, axis, and legend. Save plot.
    ax.set_xticks([p + w for p in pos], minor=False)
    ax.set_xticklabels(names, minor = False)

    plt.savefig(output)
    plt.close()
    
if __name__ == "__main__":
    main()