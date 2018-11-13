import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import common_ops as ops

"""
For each of the annotations, find its information gain for each sig.
For those sig-annotation pairs that are significant, print them to a list.
"""
BIN_SIZE = 50
def main():

    #Read in the chromHMM and annotated files.
    prediction_bed = sys.argv[1]
    prediction_or_tss_bed = sys.argv[2]
    prediction_and_tss_bed = sys.argv[3]
    ground_truth_bed = sys.argv[4]
    permuted_bed = sys.argv[5]
    plot_path = sys.argv[8]
    cell = sys.argv[7]
    training = sys.argv[6]

    #For each chromosome, get all distributions.
    predictions = []
    predictions_or = []
    predictions_and = []
    predictions_perm = []
    ground_truth = []
    chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    labels = ["Promoter", "Enhancer", "Weak", "Other", "Unknown"]
    for chr in chroms:
        predictions.append(get_distrib(prediction_bed + "/anno" + chr + ".bed"))
        predictions_perm.append(get_distrib(permuted_bed + "/anno" + chr + ".bed"))
        predictions_or.append(get_distrib(prediction_or_tss_bed + "/" + chr + ".bed"))
        predictions_and.append(get_distrib(prediction_and_tss_bed + "/anno" + chr + ".bed"))
        ground_truth.append(get_distrib(ground_truth_bed + "/chr" + chr + ".bed"))
        
    #Save reports over all chromosomes, excluding Y.
    report = open(plot_path + "_report", "w")
    report.write(",".join(str(i) for i in np.sum(np.asarray(predictions)[:22,:], axis = 0) / 22) + "\n")
    report.write(",".join(str(i) for i in np.sum(np.asarray(predictions_perm)[:22,:], axis = 0) / 22) + "\n")
    report.write(",".join(str(i) for i in np.sum(np.asarray(predictions_or)[:22,:], axis = 0) / 22) + "\n")
    report.write(",".join(str(i) for i in np.sum(np.asarray(predictions_and)[:22,:], axis = 0) / 22) + "\n")
    report.write(",".join(str(i) for i in np.sum(np.asarray(ground_truth)[:22,:], axis = 0) / 22) + "\n")
    report.close()
    
    #Plot distributions in a bar plot for all chromosomes.
    plot_bars(chroms, labels, predictions, predictions_or, predictions_and, predictions_perm, ground_truth, cell, training, plot_path)

#Get the distribution of region annotations.
def get_distrib(file):
        
    #Set up vector of total annotations.
    sum_vec = np.zeros(5)
    
    #Get scores and labels for each bed file.
    bed = np.genfromtxt(file , delimiter='\t', dtype = str)
    promoter = ["Promoter", "1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk"]
    enhancer = ["Enhancer", "6_EnhG", "7_Enh", "12_EnhBiv"]
    weak = ["Weak", "9_Het", "15_Quies"]
    other = ["Other"]
    unknown = ["Unknown"]

    #Loop through bed file and keep a running sum of each type of annotation.
    for i in range(0, bed.shape[0]):            
            
        #Get the next element data.
        next_line = bed[i,:]
        size = int(next_line[2]) - int(next_line[1])
        a = next_line[3]
       
        #Add the length of the region to the total for that region.
        if a in promoter:
            sum_vec[0] += size
        elif a in enhancer:
            sum_vec[1] += size
        elif a in weak:
            sum_vec[2] += size
        elif a in unknown:
            sum_vec[4] += size
        elif a in other or a != ".":
            sum_vec[3] += size
            
    #Compute final percentages.
    total = np.sum(sum_vec)
    total_vec = np.tile(total, 5)
    percent_vec = sum_vec / total_vec
    
    #Return value.
    return percent_vec
    
"""
Plot annotations in a stacked bar plot.
""" 
def plot_bars(names, label_names, predictions, predictions_or, predictions_and, predictions_perm, ground_truth, cell, training, output):

    #Get all percentages for each chromosome (for shape-based predictions).
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
    
    #Get all percentages for each chromosome (for tss-based OR predictions).
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
    
    #Get all percentages for each chromosome (for tss-based predictions).
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
    
    #Get all percentages for each chromosome (for permuted predictions).
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
    fig, ax = plt.subplots(figsize=(15,5))
    wid = 0.75
    w = 1
    ax.set_xlim(-1, len(pos) * bar_count)
    color_unknown = "black"
    color_other = "gray"
    color_promoter = "green"
    color_enhancer = "blue"
    color_weak = "red"
    
    #Plot shape-based predictions.
    p1_s = plt.bar(pos, unknowns_s, color = color_unknown, width = wid)
    p2_s = plt.bar(pos, other_s, color = color_other, bottom = unknowns_s, width = wid)
    p3_s = plt.bar(pos, promoters_s, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_s, other_s)],  width = wid)
    p4_s = plt.bar(pos, enhancers_s, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_s, other_s, promoters_s)], width = wid)
    p5_s = plt.bar(pos, weaks_s, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_s, other_s, promoters_s, enhancers_s)], width = wid)
    for p in pos:
        ax.text(p, 1.0, "S", ha='center', va='bottom')
    
    #Plot shape-based predictions OR position-based predictions.
    pos2 = [p + w for p in pos]
    p1_t = plt.bar(pos2, unknowns_or, color = color_unknown, width = wid)
    p2_t = plt.bar(pos2, other_or, color = color_other, bottom = unknowns_or, width = wid)
    p3_t = plt.bar(pos2, promoters_or, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_or, other_or)],  width = wid)
    p4_t = plt.bar(pos2, enhancers_or, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_or, other_or, promoters_or)], width = wid)
    p5_t = plt.bar(pos2, weaks_or, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_or, other_or, promoters_or, enhancers_or)], width = wid)
    for p in pos2:
        ax.text(p, 1.0, "X", ha='center', va='bottom')
        
    #Plot shape-based predictions OR position-based predictions.
    pos3 = [p + w * 2 for p in pos]
    p1_ta = plt.bar(pos3, unknowns_and, color = color_unknown, width = wid)
    p2_ta = plt.bar(pos3, other_and, color = color_other, bottom = unknowns_and, width = wid)
    p3_ta = plt.bar(pos3, promoters_and, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_and, other_and)],  width = wid)
    p4_ta = plt.bar(pos3, enhancers_and, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_and, other_and, promoters_and)], width = wid)
    p5_ta = plt.bar(pos3, weaks_and, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_and, other_and, promoters_and, enhancers_and)], width = wid)
    for p in pos3:
        ax.text(p, 1.0, "+", ha='center', va='bottom')
    
    #Plot ground truth.
    pos4 = [p + w * 3 for p in pos]
    p1_g = plt.bar(pos4, unknowns_g, color = color_unknown, width = wid, label = "Unknown")
    p2_g = plt.bar(pos4, other_g, color = color_other, bottom = unknowns_g, width = wid, label = "Other")
    p3_g = plt.bar(pos4, promoters_g, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_g, other_g)],  width = wid, label = "Promoter")
    p4_g = plt.bar(pos4, enhancers_g, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_g, other_g, promoters_g)], width = wid, label = "Enhancer")
    p5_g = plt.bar(pos4, weaks_g, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_g, other_g, promoters_g, enhancers_g)], width = wid, label = "Weak")
    for p in pos4:
        ax.text(p, 1.0, "G", ha='center', va='bottom')
        
    #Plot permuted predictions.
    # pos4 = [p + w * 3 for p in pos]
    # p1_p = plt.bar(pos4, unknowns_p, color = color_unknown, width = wid)
    # p2_p = plt.bar(pos4, other_p, color = color_other, bottom = unknowns_p, width = wid)
    # p3_p = plt.bar(pos4, promoters_p, color = color_promoter, bottom = [i+j for i,j in zip(unknowns_p, other_p)],  width = wid)
    # p4_p = plt.bar(pos4, enhancers_p, color = color_enhancer, bottom = [i+j+k for i,j,k in zip(unknowns_p, other_p, promoters_p)], width = wid)
    # p5_p = plt.bar(pos4, weaks_p, color = color_weak, bottom = [i+j+k+l for i,j,k,l in zip(unknowns_p, other_p, promoters_p, enhancers_p)], width = wid)
    # for p in pos4:
        # ax.text(p, 1.0, "P", ha='center', va='bottom')
    
    plt.title("Annotation Distribution for All Chromosomes - " + cell + " from " + training)
    plt.xlabel("Chromosome")
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
                        Patch(facecolor=color_unknown, label='Unknown'),
                        Patch(facecolor=color_other, label='Other')
                       ]
    plt.legend(handles=legend_elements, loc="lower right")
    plt.savefig(output + "_legend" + ".png")

    
if __name__ == "__main__":
    main()