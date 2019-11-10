#Import sys for obtaining command line args.
import numpy as np
import pandas as pd

#Import matplotlib for plots.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn

# Constants
CELL_TYPES = ["GM12878", "A549", "Brain", "H1"]
CHROMOSOMES = [str(i) for i in range(1,23)]
REPEATS = [str(i) for i in range(1,6)]
ALPHAS = ['{0:3.1f}'.format(i/10.0) for i in range(1,10)]
SIGMAS = [str(i) for i in range(1,10)]
K = ["20", "40", "80", "160"]
MAX_DIST = ['{0:3.1f}'.format(i/10.0) for i in range(1,10)]

"""
The following functions are for consolidating all precision and recall methods for a given
method and experiment.
"""
def precision_recall_som_intrachrom(directory, avg_over_alpha, avg_over_sigma, avg_over_repeats, avg_over_cell_type, avg_over_chromosome):

    # Dimensions: cell types, chromosomes, repeats, alpha, sigma, RE, precision/recall
    all_precision_recall = np.zeros((4,22,5,9,9,4,2))
    
    for i in range(len(CELL_TYPES)):
        for j in range(len(CHROMOSOMES)):
            for k in range(len(REPEATS)):
                for l in range(len(ALPHAS)):
                    for m in range(len(SIGMAS)):
                    
                        precision_recall = pd.read_csv(open(directory + "/" + CELL_TYPES[i] + "/" + CHROMOSOMES[i] + "/" + REPEATS[i] + "/" + ALPHAS[i] + "/" + SIGMAS[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
                        all_precision_recall[i,j,k,l,m,:,:] = np.asarray(precision_recall)
                        
    # Average over requested parameters.
    if avg_over_sigma:
        all_precision_recall = np.mean(all_precision_recall, axis = 4)
    if avg_over_alpha:
        all_precision_recall = np.mean(all_precision_recall, axis = 3)
    if avg_over_repeats:
        all_precision_recall = np.mean(all_precision_recall, axis = 2)
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
                        
    
    # $BASE_PATH/som_vn_intrachrom/
    # $BASE_PATH/som_vn_peas_intrachrom/
    # $BASE_PATH/som_intrachrom/
    # $BASE_PATH/som_peas_intrachrom/
    # $BASE_PATH/som_vn_signalperm_intrachrom/
    # $BASE_PATH/som_vn_signalperm_peas_intrachrom/
    # $BASE_PATH/som_vn_chromhmmperm_intrachrom/
    # $BASE_PATH/som_vn_chromhmmperm_peas_intrachrom/
    # ${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
    
def precision_recall_cagt_intrachrom(directory, avg_over_max_dist, avg_over_k, avg_over_repeats, avg_over_chromosome, avg_over_cell_type):
    
    # Dimensions: cell types, chromosomes, repeats, k, max_dist, RE, precision/recall
    all_precision_recall = np.zeros((4,22,5,4,9,4,2))
    
    for i in range(len(CELL_TYPES)):
        for j in range(len(CHROMOSOMES)):
            for k in range(len(REPEATS)):
                for l in range(len(K)):
                    for m in range(len(MAX_DIST)):
                    
                        precision_recall = pd.read_csv(open(directory + "/" + CELL_TYPES[i] + "/" + CHROMOSOMES[i] + "/" + REPEATS[i] + "/" + K[i] + "/" + MAX_DIST[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
                        all_precision_recall[i,j,k,l,m,:,:] = np.asarray(precision_recall)
                        
    # Average over requested parameters.
    if avg_over_k:
        all_precision_recall = np.mean(all_precision_recall, axis = 4)
    if avg_over_max_dist:
        all_precision_recall = np.mean(all_precision_recall, axis = 3)
    if avg_over_repeats:
        all_precision_recall = np.mean(all_precision_recall, axis = 2)
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/cagt_intrachrom/
    # $BASE_PATH/cagt_peas_intrachrom/
    # ${cell_types[$i]}/$chrom/$repeat/$k/$max_dist/
    
def precision_recall_signal_intrachrom(directory, avg_over_chromosome, avg_over_cell_type):
    
    # Dimensions: cell types, chromosomes, RE, precision/recall
    all_precision_recall = np.zeros((4,22,4,2))
    
    for i in range(len(CELL_TYPES)):
        for j in range(len(CHROMOSOMES)):
                    
            precision_recall = pd.read_csv(open(directory + "/" + CELL_TYPES[i] + "/" + CHROMOSOMES[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
            all_precision_recall[i,j,:,:] = np.asarray(precision_recall)
            
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/signal_intrachrom/
    # $BASE_PATH/signal_peas_intrachrom/
    # ${cell_types[$i]}/$chrom/
    
def precision_recall_som_crosschrom(directory, avg_over_alpha, avg_over_sigma, avg_over_repeats, avg_over_chromosome, avg_over_cell_type):
    
    # Dimensions: cell types, repeats, alpha, sigma, RE, precision/recall
    all_precision_recall = np.zeros((4,5,9,9,4,2))
    
    for i in range(len(CELL_TYPES)):
        for k in range(len(REPEATS)):
            for l in range(len(ALPHAS)):
                for m in range(len(SIGMAS)):
                
                    precision_recall = pd.read_csv(open(directory + "/" + CELL_TYPES[i] + "/" + REPEATS[i] + "/" + ALPHAS[i] + "/" + SIGMAS[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
                    all_precision_recall[i,k,l,m,:,:] = np.asarray(precision_recall)
                    
    # Average over requested parameters.
    if avg_over_sigma:
        all_precision_recall = np.mean(all_precision_recall, axis = 4)
    if avg_over_alpha:
        all_precision_recall = np.mean(all_precision_recall, axis = 3)
    if avg_over_repeats:
        all_precision_recall = np.mean(all_precision_recall, axis = 2)
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/som_vn_${chrom_i}_${chrom_j}/
    # $BASE_PATH/som_vn_${chrom_i}_${chrom_j}_peas/
    # $BASE_PATH/som_${chrom_i}_${chrom_j}/
    # $BASE_PATH/som_${chrom_i}_${chrom_j}_peas/
    # $BASE_PATH/som_vn_signalperm_${chrom_i}_${chrom_j}/
    # $BASE_PATH/som_vn_signalperm_${chrom_i}_${chrom_j}_peas/
    # $BASE_PATH/som_vn_chromhmmperm_${chrom_i}_${chrom_j}/
    # $BASE_PATH/som_vn_chromhmmperm_${chrom_i}_${chrom_j}_peas/
    # ${cell_types[$i]}/$repeat/$alpha_0/$sigma_0/
    
def precision_recall_cagt_crosschrom(directory, avg_over_max_dist, avg_over_k, avg_over_repeats):
    
    # Dimensions: cell types, repeats, k, max_dist, RE, precision/recall
    all_precision_recall = np.zeros((4,5,4,9,4,2))
    
    for i in range(len(CELL_TYPES)):
        for k in range(len(REPEATS)):
            for l in range(len(K)):
                for m in range(len(MAX_DIST)):
                
                    precision_recall = precision_recall = pd.read_csv(open(directory + "/" + CELL_TYPES[i] + "/" + REPEATS[i] + "/" + K[i] + "/" + MAX_DIST[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
                    all_precision_recall[i,k,l,m,:,:] = np.asarray(precision_recall)
                    
    # Average over requested parameters.
    if avg_over_k:
        all_precision_recall = np.mean(all_precision_recall, axis = 4)
    if avg_over_max_dist:
        all_precision_recall = np.mean(all_precision_recall, axis = 3)
    if avg_over_repeats:
        all_precision_recall = np.mean(all_precision_recall, axis = 2)
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/cagt_intrachrom/
    # $BASE_PATH/cagt_peas_intrachrom/
    # ${cell_types[$i]}/$chrom/$repeat/$k/$max_dist/
    
def precision_recall_signal_crosschrom(directory, avg_over_chromosome, avg_over_cell_type):
    
    # Dimensions: cell types, RE, precision/recall
    all_precision_recall = np.zeros((4,4,2))
    
    for i in range(len(CELL_TYPES)):
        for j in range(len(CHROMOSOMES)):
                    
            precision_recall = precision_recall = pd.read_csv(open(directory + "/" + CELL_TYPES[i] + "/" + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
            all_precision_recall[i,j,:,:] = np.asarray(precision_recall)
            
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/signal_intrachrom/
    # $BASE_PATH/signal_peas_intrachrom/
    # ${cell_types[$i]}/$chrom/
    
def precision_recall_som_crosscelltype(directory, avg_over_alpha, avg_over_sigma, avg_over_repeats, avg_over_cell_type, avg_over_chromosome):

    # Dimensions: cell types, chromosomes, repeats, alpha, sigma, RE, precision/recall
    all_precision_recall = np.zeros((22,5,9,9,4,2))
    
    for j in range(len(CHROMOSOMES)):
        for k in range(len(REPEATS)):
            for l in range(len(ALPHAS)):
                for m in range(len(SIGMAS)):
                
                    precision_recall = precision_recall = pd.read_csv(open(directory + "/" + CHROMOSOMES[i] + "/" + REPEATS[i] + "/" + ALPHAS[i] + "/" + SIGMAS[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
                    all_precision_recall[j,k,l,m,:,:] = np.asarray(precision_recall)
                    
    # Average over requested parameters.
    if avg_over_sigma:
        all_precision_recall = np.mean(all_precision_recall, axis = 4)
    if avg_over_alpha:
        all_precision_recall = np.mean(all_precision_recall, axis = 3)
    if avg_over_repeats:
        all_precision_recall = np.mean(all_precision_recall, axis = 2)
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/som_vn_${cell_types[$i]}_${cell_types[$j]}/
    # $BASE_PATH/som_vn_${cell_types[$i]}_${cell_types[$j]}_peas/
    # $BASE_PATH/som_${cell_types[$i]}_${cell_types[$j]}/
    # $BASE_PATH/som_${cell_types[$i]}_${cell_types[$j]}_peas/
    # $BASE_PATH/som_vn_${cell_types[$i]}_${cell_types[$j]}_signalperm/
    # $BASE_PATH/som_vn_${cell_types[$i]}_${cell_types[$j]}_signalperm_peas/
    # $BASE_PATH/som_vn_${cell_types[$i]}_${cell_types[$j]}_chromhmmperm/
    # $BASE_PATH/som_vn_${cell_types[$i]}_${cell_types[$j]}_chromhmmperm_peas/
    # $chrom/$repeat/$alpha_0/$sigma_0/
    
def precision_recall_cagt_interchrom(directory, avg_over_max_dist, avg_over_k, avg_over_repeats, avg_over_chromosome):
    
    # Dimensions: cell types, chromosomes, repeats, k, max_dist, RE, precision/recall
    all_precision_recall = np.zeros((22,5,4,9,4,2))
    
    for j in range(len(CHROMOSOMES)):
        for k in range(len(REPEATS)):
            for l in range(len(K)):
                for m in range(len(MAX_DIST)):
                
                    precision_recall = precision_recall = pd.read_csv(open(directory + "/" + CHROMOSOMES[i] + "/" + REPEATS[i] + "/" + K[i] + "/" + MAX_DIST[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
                    all_precision_recall[j,k,l,m,:,:] = np.asarray(precision_recall)
                    
    # Average over requested parameters.
    if avg_over_k:
        all_precision_recall = np.mean(all_precision_recall, axis = 4)
    if avg_over_max_dist:
        all_precision_recall = np.mean(all_precision_recall, axis = 3)
    if avg_over_repeats:
        all_precision_recall = np.mean(all_precision_recall, axis = 2)
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/cagt_${cell_types[$i]}_${cell_types[$j]}/
    # $BASE_PATH/cagt_${cell_types[$i]}_${cell_types[$j]}_peas/
    # $chrom/$repeat/$k/$max_dist/
    
def precision_recall_signal_interchrom(directory, avg_over_cell_type, avg_over_chromosome):
    
    # Dimensions: cell types, chromosomes, RE, precision/recall
    all_precision_recall = np.zeros((22,4,2))
    
    for j in range(len(CHROMOSOMES)):
                
        precision_recall = precision_recall = pd.read_csv(open(directory + "/" + CHROMOSOMES[i] + "/" + REPEATS[i] + "/" + ALPHAS[i] + "/" + SIGMAS[i] + "/precision_recall" + CHROMOSOMES[i] + ".csv", delim = "\t")
        all_precision_recall[j,:,:] = np.asarray(precision_recall)
        
    if avg_over_chromosome:
        all_precision_recall = np.mean(all_precision_recall, axis = 1)
    if avg_over_cell_type:
        all_precision_recall = np.mean(all_precision_recall, axis = 0)
        
    return all_precision_recall
    
    # $BASE_PATH/signal_${cell_types[$i]}_${cell_types[$j]}/
    # $BASE_PATH/signal_${cell_types[$i]}_${cell_types[$j]}_peas/
    # $chrom/

"""
Plot the percentage distribution for
each type of RE.
"""
def make_one_plot(title, som_vn, som, cagt, signal, som_vn_signalperm, som_vn_chromhmmperm):

    # Get AUC for all models.
    som_vn_auc = sklearn.metrics.auc(som_vn[0], som_vn[1])
    som_auc = sklearn.metrics.auc(som[0], som[1])
    cagt_auc = sklearn.metrics.auc(cagt[0], cagt[1])
    signal_auc = sklearn.metrics.auc(signal[0], signal[1])
    som_vn_signalperm_auc = sklearn.metrics.auc(som_vn_signalperm[0], som_vn_signalperm[1])
    som_vn_chromhmmperm_auc = sklearn.metrics.auc(som_vn_chromhmmperm[0], som_vn_chromhmmperm[1])

    # Plot each model's curve.
    plt.plot(som_vn[0], som_vn[1], label='SOM-VN (area = " + str(som_vn_auc))
    plt.plot(som[0], som[1]], label='SOM (area = " + str(som_auc))
    plt.plot(cagt[0], cagt[1]], label='CAGT (area = " + str(cagt_auc))
    plt.plot(signal[0], signal[1]], label='Signal Intensity (area = " + str(signal_auc))
    plt.plot(som_vn_signalperm[0], som_vn_signalperm[1], label='SOM-VN (area = " + str(som_vn_signalperm_auc))
    plt.plot(som_vn_chromhmmperm[0], som_vn_chromhmmperm[1], label='SOM-VN (area = " + str(som_vn_chromhmmperm_auc))

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title)
    plt.legend(loc="lower right")

if __name__ == "__main__":
    main()