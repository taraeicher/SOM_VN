import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({'font.size': 14})

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #For each chromosome, get all distributions.
    files = []
    combo = []
    cells = ["A549", "Brain", "H1"]
    for src in cells:
        for dest in cells:
            files.append(sys.argv[1] + dest + sys.argv[2] + "/" + src + "_to_" + dest + ".png_report")
            combo.append(src[0] + "-" + dest[0])
    plot_out = sys.argv[3]
    
    #Get all distributions.
    predictions = []
    predictions_or = []
    predictions_and = []
    predictions_perm = []
    ground_truth = []
    labels = ["Promoter", "Enhancer", "Weak", "Other", "Unknown"]
    for f in files:
        file = open(f, "r")
        predictions.append([float(s) for s in file.readline().split(",")])
        predictions_perm.append([float(s) for s in file.readline().split(",")])
        predictions_or.append([float(s) for s in file.readline().split(",")])
        predictions_and.append([float(s) for s in file.readline().split(",")])
        ground_truth.append([float(s) for s in file.readline().split(",")])
    
    #Plot distributions in a bar plot for all cell combos.
    plot_bars(combo, labels, predictions, predictions_or, predictions_and, predictions_perm, ground_truth, plot_out)
    
"""
Plot annotations in a stacked bar plot.
""" 
def plot_bars(names, label_names, predictions, predictions_or, predictions_and, predictions_perm, ground_truth, output):

    #Get all percentages for each combination (for shape-based predictions).
    unknown_percentage_s = []
    promoter_percentage_s = []
    enhancer_percentage_s = []
    weak_percentage_s = []
    other_percentage_s = []
    for c in range(0, len(names)):
        unknown_percentage_s.append(predictions[c][label_names.index("Unknown")])
        promoter_percentage_s.append(predictions[c][label_names.index("Promoter")])
        enhancer_percentage_s.append(predictions[c][label_names.index("Enhancer")])
        weak_percentage_s.append(predictions[c][label_names.index("Weak")])
        other_percentage_s.append(predictions[c][label_names.index("Other")])
    unknowns_s = np.asarray([float(i) for i in unknown_percentage_s])
    promoters_s = np.asarray([float(i) for i in promoter_percentage_s])
    enhancers_s = np.asarray([float(i) for i in enhancer_percentage_s])
    weaks_s = np.asarray([float(i) for i in weak_percentage_s])
    other_s = np.asarray([float(i) for i in other_percentage_s])
    
    #Get all percentages for each combination (for tss-based OR predictions).
    unknown_percentage_or = []
    promoter_percentage_or = []
    enhancer_percentage_or = []
    weak_percentage_or = []
    other_percentage_or = []
    for c in range(0, len(names)):
        unknown_percentage_or.append(predictions_or[c][label_names.index("Unknown")])
        promoter_percentage_or.append(predictions_or[c][label_names.index("Promoter")])
        enhancer_percentage_or.append(predictions_or[c][label_names.index("Enhancer")])
        weak_percentage_or.append(predictions_or[c][label_names.index("Weak")])
        other_percentage_or.append(predictions_or[c][label_names.index("Other")])
    unknowns_or = np.asarray([float(i) for i in unknown_percentage_or])
    promoters_or = np.asarray([float(i) for i in promoter_percentage_or])
    enhancers_or = np.asarray([float(i) for i in enhancer_percentage_or])
    weaks_or = np.asarray([float(i) for i in weak_percentage_or])
    other_or = np.asarray([float(i) for i in other_percentage_or])
    
    #Get all percentages for each combination (for tss-based AND predictions).
    unknown_percentage_and = []
    promoter_percentage_and = []
    enhancer_percentage_and = []
    weak_percentage_and = []
    other_percentage_and = []
    for c in range(0, len(names)):
        unknown_percentage_and.append(predictions_and[c][label_names.index("Unknown")])
        promoter_percentage_and.append(predictions_and[c][label_names.index("Promoter")])
        enhancer_percentage_and.append(predictions_and[c][label_names.index("Enhancer")])
        weak_percentage_and.append(predictions_and[c][label_names.index("Weak")])
        other_percentage_and.append(predictions_and[c][label_names.index("Other")])
    unknowns_and = np.asarray([float(i) for i in unknown_percentage_and])
    promoters_and = np.asarray([float(i) for i in promoter_percentage_and])
    enhancers_and = np.asarray([float(i) for i in enhancer_percentage_and])
    weaks_and = np.asarray([float(i) for i in weak_percentage_and])
    other_and = np.asarray([float(i) for i in other_percentage_and])
    
    #Get all percentages for each combination (for permuted predictions).
    unknown_percentage_p = []
    promoter_percentage_p = []
    enhancer_percentage_p = []
    weak_percentage_p = []
    other_percentage_p = []
    for c in range(0, len(names)):
        unknown_percentage_p.append(predictions_perm[c][label_names.index("Unknown")])
        promoter_percentage_p.append(predictions_perm[c][label_names.index("Promoter")])
        enhancer_percentage_p.append(predictions_perm[c][label_names.index("Enhancer")])
        weak_percentage_p.append(predictions_perm[c][label_names.index("Weak")])
        other_percentage_p.append(predictions_perm[c][label_names.index("Other")])
    unknowns_p = np.asarray([float(i) for i in unknown_percentage_p])
    promoters_p = np.asarray([float(i) for i in promoter_percentage_p])
    enhancers_p = np.asarray([float(i) for i in enhancer_percentage_p])
    weaks_p = np.asarray([float(i) for i in weak_percentage_p])
    other_p = np.asarray([float(i) for i in other_percentage_p])
    
    #Get all percentages for each chromosome (for ChromHMM ground truth).
    unknown_percentage_g = []
    promoter_percentage_g = []
    enhancer_percentage_g = []
    weak_percentage_g = []
    other_percentage_g = []
    for c in range(0, len(names)):
        unknown_percentage_g.append(ground_truth[c][label_names.index("Unknown")])
        promoter_percentage_g.append(ground_truth[c][label_names.index("Promoter")])
        enhancer_percentage_g.append(ground_truth[c][label_names.index("Enhancer")])
        weak_percentage_g.append(ground_truth[c][label_names.index("Weak")])
        other_percentage_g.append(ground_truth[c][label_names.index("Other")])
    unknowns_g = np.asarray([float(i) for i in unknown_percentage_g])
    promoters_g = np.asarray([float(i) for i in promoter_percentage_g])
    enhancers_g = np.asarray([float(i) for i in enhancer_percentage_g])
    weaks_g = np.asarray([float(i) for i in weak_percentage_g])
    other_g = np.asarray([float(i) for i in other_percentage_g])

    #Set up parameters for plot.
    bar_count = 5
    pos = list(range(0, len(unknowns_g) * bar_count, bar_count)) 
    fig, ax = plt.subplots(figsize=(10,5))
    wid = 0.75
    w = 1
    ax.set_xlim(-1, len(pos) * bar_count)
    color_unknown = "white"
    color_other = "white"
    color_promoter = "black"
    color_enhancer = "gray"
    color_weak = "silver"
    
    #Plot shape-based predictions.
    p1_s = plt.bar(pos, unknowns_s, color = color_unknown, edgecolor = "black", width = wid)
    p2_s = plt.bar(pos, other_s, color = color_other, bottom = unknowns_s, width = wid)
    p3_s = plt.bar(pos, promoters_s, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_s, other_s)],  width = wid)
    p4_s = plt.bar(pos, enhancers_s, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_s, other_s, promoters_s)], width = wid)
    p5_s = plt.bar(pos, weaks_s, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_s, other_s, promoters_s, enhancers_s)], width = wid)
    for p in pos:
        ax.text(p, 1.0, "S", ha='center', va='bottom')
    
    #Plot shape-based predictions OR position-based predictions.
    pos2 = [p + w for p in pos]
    p1_t = plt.bar(pos2, unknowns_or, color = color_unknown, edgecolor = "black", width = wid)
    p2_t = plt.bar(pos2, other_or, color = color_other, edgecolor = "black", bottom = unknowns_or, width = wid)
    p3_t = plt.bar(pos2, promoters_or, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_or, other_or)],  width = wid)
    p4_t = plt.bar(pos2, enhancers_or, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_or, other_or, promoters_or)], width = wid)
    p5_t = plt.bar(pos2, weaks_or, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_or, other_or, promoters_or, enhancers_or)], width = wid)
    for p in pos2:
        ax.text(p, 1.0, "X", ha='center', va='bottom')
        
    #Plot shape-based predictions AND position-based predictions.
    pos3 = [p + w * 2 for p in pos]
    p1_t = plt.bar(pos3, unknowns_and, color = color_unknown, edgecolor = "black", width = wid)
    p2_t = plt.bar(pos3, other_and, color = color_other, edgecolor = "black", bottom = unknowns_and, width = wid)
    p3_t = plt.bar(pos3, promoters_and, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_and, other_and)],  width = wid)
    p4_t = plt.bar(pos3, enhancers_and, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_and, other_and, promoters_and)], width = wid)
    p5_t = plt.bar(pos3, weaks_and, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_and, other_and, promoters_and, enhancers_and)], width = wid)
    for p in pos3:
        ax.text(p, 1.0, "+", ha='center', va='bottom')
    
    #Plot ground truth.
    pos4 = [p + w * 3 for p in pos]
    p1_g = plt.bar(pos4, unknowns_g, color = color_unknown, edgecolor = "black", width = wid, label = "Unknown")
    p2_g = plt.bar(pos4, other_g, color = color_other, edgecolor = "black", bottom = unknowns_g, width = wid, label = "Other")
    p3_g = plt.bar(pos4, promoters_g, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_g, other_g)],  width = wid, label = "Promoter")
    p4_g = plt.bar(pos4, enhancers_g, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_g, other_g, promoters_g)], width = wid, label = "Enhancer")
    p5_g = plt.bar(pos4, weaks_g, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_g, other_g, promoters_g, enhancers_g)], width = wid, label = "Weak")
    for p in pos4:
        ax.text(p, 1.0, "G", ha='center', va='bottom')

    plt.xlabel("Training and Testing Pair")
    plt.ylabel("Percentage") 

    # Add title, axis, and legend. Save plot.
    ax.set_xticks([p + w for p in pos], minor=False)
    ax.set_xticklabels(names, minor = False)
    
    plt.savefig(output + ".png")
    plt.close()
    
    #Save the legend.
    legend_elements = [Patch(facecolor=color_enhancer, label='Enhancer'),
                        Patch(facecolor=color_promoter, label='Promoter'),
                        Patch(facecolor=color_weak, label='Weak'),
                        Patch(facecolor=color_unknown, edgecolor = "black", label='Unknown / Other'),
                       ]
    plt.legend(handles=legend_elements, loc="lower right")
    plt.savefig(output + "_legend" + ".png")
    
if __name__ == "__main__":
    main()