from sklearn.cluster import KMeans
import numpy as np
import scipy as sp
import sys
import os
import imp
import math
import common_ops as ops

gap = imp.load_source('gap', '/users/PAS0272/osu5316/miniconda3/lib/python3.5/site-packages/gap/gap.py')
#gap_src.gap()
#from gap import gap

"""
For each of the annotations, find its information gain for each cluster.
For those cluster-annotation pairs that are significant, print them to a list.
"""
def main():

	#Read in the bed file and cluster data.
	bed = np.genfromtxt(sys.argv[1], delimiter='\t', dtype = str)
	#bed = bed[np.where(bed[:,4] != '.')]
	clusters = []
	file = open(sys.argv[2], 'r')
	output = sys.argv[3]
	file_annotation = sys.argv[4]
	wig_file = sys.argv[5]
	next_line = file.readline()
	while next_line:
		clusters.append(next_line)
		next_line = file.readline()
		
	#Cluster annotations are in column 3. Biological annotations are in cluster 4.
	cluster_col = 3
	bio_col = 4#8
	
	#Merge repeated regions together.
	merge_list(bed[:,bio_col])
	
	#Find the list of unique clusters and annotations.
	unique_clusts = list(set(bed[:,cluster_col]))
	unique_anno = list(set(bed[:,bio_col]))

	#Using the cardinality and entropy values for all annotations and cluster-annotation pairs,
	#calculate information gain for each cluster with respect to each annotation.
	anno_probs = get_cardinalities(bed, unique_anno, bio_col) / len(bed)
	cluster_probs = get_cluster_cardinalities(bed, unique_clusts, cluster_col) / len(bed)
	all_cards, all_op_anno_cards, all_op_clust_cards, all_op_cards = get_all_cardinalities(bed, unique_anno, unique_clusts, cluster_col, bio_col)
	info_gains, in_cluster_sig = get_adj_info_gains(anno_probs, cluster_probs, all_cards / len(bed), all_op_anno_cards / len(bed), all_op_clust_cards / len(bed), all_op_cards / len(bed))
	
	#Print all clusters with significant annotations, along with their annotations.
	print_gains_above_threshold(info_gains, 0.6, clusters, unique_clusts, unique_anno, output, file_annotation, in_cluster_sig, wig_file)

#Each row in the bio column is expected to contain a list of annotations present in the cluster region.
#This method removes repeats in the list.
def merge_list(lst_vector):
	for i in range(0, len(lst_vector)):
		terms = lst_vector[i].split(',')
		
		#Replace terms with their equivalent biological classes.
		j = 0
		while j < len(terms):
			if terms[j] == "1_TssA" or terms[j] == "2_TssAFlnk" or terms[j] == "10_TssBiv" or terms[j] == "11_BivFlnk":
				terms[j] = "Promoter"
				j = j + 1
			elif terms[j] == "6_EnhG" or terms[j] == "7_Enh" or terms[j] == "12_EnhBiv":
				terms[j] = "Enhancer"
				j = j + 1
			elif terms[j] == "13_ReprPC" or terms[j] == "ReprPCWk":
				terms[j] = "Polycomb"
				j = j + 1
			elif terms[j] == "9_Het" or terms[j] == "15_Quies":
				terms[j] = "Weak"
				j = j + 1
			else:
				terms.pop(j)
				
		#Joins elements.
		lst_vector[i] = ','.join(sorted(list(set(terms))))
	
#Find the cardinalities of all unique annotations.
def get_cardinalities(bed, anno_list, col):
	cardinalities = np.zeros(len(anno_list))
	for i in range(0, len(anno_list)):
		cardinalities[i] = len(np.where(bed[:,col] == anno_list[i])[0])
	return cardinalities

#Find the cardinalities of all unique clusters.
def get_cluster_cardinalities(bed, clust_list, col):
	cardinalities = np.zeros(len(clust_list))
	for i in range(0, len(clust_list)):
		cardinalities[i] = len(np.where(bed[:,col] == clust_list[i])[0])
	return cardinalities
	
#Find the cardinalities of all unique clusters given the annotation.
def get_all_cardinalities(bed, anno_list, clust_list, col, bio_col):
	cardinalities = np.zeros((len(clust_list), len(anno_list)))
	opp_anno_cardinalities = np.zeros((len(clust_list), len(anno_list)))
	opp_clust_cardinalities = np.zeros((len(clust_list), len(anno_list)))
	opp_both_cardinalities = np.zeros((len(clust_list), len(anno_list)))
	for i in range(0, len(clust_list)):
		for j in range(0, len(anno_list)):
			cardinalities[i, j] = len(bed[list(set(list(np.where(bed[:,col] == clust_list[i])[0])) & set(list(np.where(bed[:, bio_col] == anno_list[j])[0]))),:])
			opp_anno_cardinalities[i, j] = len(bed[list(set(list(np.where(bed[:,col] == clust_list[i])[0])) & set(list(np.where(bed[:, bio_col] != anno_list[j])[0]))),:])
			opp_clust_cardinalities[i, j] = len(bed[list(set(list(np.where(bed[:,col] != clust_list[i])[0])) & set(list(np.where(bed[:, bio_col] == anno_list[j])[0]))),:])
			opp_both_cardinalities[i, j] = len(bed[list(set(list(np.where(bed[:,col] != clust_list[i])[0])) & set(list(np.where(bed[:, bio_col] != anno_list[j])[0]))),:])
	return cardinalities, opp_anno_cardinalities, opp_clust_cardinalities, opp_both_cardinalities

#Find the info gains adjusted by annotation probabilities.
def get_adj_info_gains(anno_probs, cluster_probs, all_probs, all_op_anno_probs, all_op_clust_probs, all_op_probs):
	probs_all = np.reshape(np.tile(1 / anno_probs, len(cluster_probs)), (len(cluster_probs), len(anno_probs)))
	base_ig, in_cluster_assoc = get_info_gains(anno_probs, cluster_probs, all_probs, all_op_anno_probs, all_op_clust_probs, all_op_probs)
	return probs_all * base_ig, in_cluster_assoc

#Find the info gains adjusted by annotation probabilities.
def get_info_gains(anno_probs, cluster_probs, all_probs, all_op_anno_probs, all_op_clust_probs, all_op_probs):

	#Calculate expected cluster-annotation probability.
	cross_probs = np.transpose(np.outer(anno_probs, cluster_probs))
	opp_anno_cross_probs = np.transpose(np.outer(1 - anno_probs, cluster_probs))
	opp_clust_cross_probs = np.transpose(np.outer(anno_probs, 1 - cluster_probs))
	opp_cross_probs = np.transpose(np.outer(1 - anno_probs, 1 - cluster_probs))
	
	#Calculate ratio of actual cluster-annotation probability to expected cluster-annotation probability.
	prob_ratios = all_probs / cross_probs
	prob_ratios_op_anno = all_op_anno_probs / opp_anno_cross_probs
	prob_ratios_op_clust = all_op_clust_probs / opp_clust_cross_probs
	prob_ratios_op_all = all_op_probs / opp_cross_probs
	
	#Calculate logarithmic terms of ratios.
	log_probs = all_probs * np.log2(prob_ratios)
	log_probs_op_anno = all_op_anno_probs * np.log2(prob_ratios_op_anno)
	log_probs_op_clust = all_op_clust_probs * np.log2(prob_ratios_op_clust)
	log_probs_op = all_op_probs * np.log2(prob_ratios_op_all)
	
	#Sum together logarithmic terms for:
	#1. Has cluster annotation t and biological annotation c.
	#2. Doesn't have cluster annotation t but has biological annotation c.
	#3. Has cluster annotation t but not biological annotation c.
	#4. Doesn't have cluster annotation t or biological annotation c.
	ig = np.nan_to_num(np.add(np.add(log_probs, log_probs_op_clust), np.add(log_probs_op_anno, log_probs_op)))
	
	#Determine which is more significant: cluster t with annotation c or not cluster t with annotation c.
	anno_true = sp.greater(log_probs, log_probs_op_clust)
	
	#Return the information gains and the indicators of significance.
	return ig, anno_true
	
#Print out the cluster centroids
def print_all(info_gains, top_3, anno, file, file_annotation):
	print("\t".join(top_3[0,:]))
	file = open(file, "a")
	file.write(file_annotation + "\n")
	file.write(','.join(anno))
	file.write("\n")
	file.write(','.join(map(str, info_gains)))
	file.write("\n")
	file.write("\t".join(top_3[0,:]))
	file.write("\n")
	file.write("\t".join(top_3[1,:]))
	file.write("\n")
	file.write("\t".join(top_3[2,:]))
	file.write("\n")
	file.close()

def print_gains_above_threshold(info_gains, threshold, clusters, unique_clusts, unique_anno, filename, file_anno, in_cluster_sig, wig_file):
	
	#Open the output file.
	file = open(filename, 'a')
	wfile = open(wig_file, 'r')
	
	#Get the intensity to scale.
	intensity = ops.get_intensity_percentile(0.995, wfile)
	wfile.close()
	scale = 5.0 / intensity
	
	where_greater = np.where(np.logical_and(info_gains > threshold, in_cluster_sig))
	old_clust = "-1"
	anno_sets = list(set())
	for i in range(0,len(where_greater[0])):
		#Get the cluster and annotation split pair.
		new_clust = unique_clusts[where_greater[0][i]]
		split_anno = unique_anno[where_greater[1][i]].split(",")

		#If this is the same cluster as the previous line, keep merging annotations.
		if new_clust == old_clust:
			#If there are no elements in common anymore, start a new merged annotation.
			#Otherwise, merge again.
			j = 0
			matched = False
			while j < len(anno_sets):
				new_set = anno_sets[j].intersection(set(split_anno))
				if len(new_set) != 0:
					anno_sets[j] = new_set
					matched = True
					j = len(anno_sets)
				else:
					j = j + 1
			if matched == False:
				anno_sets.append(set(split_anno))
		
		#If this is a new cluster, start a new set of merged annotations.
		else:
			#If this is not the first cluster, print the old one with all old merged annotations.
			if old_clust != "-1":
				for j in range(0, len(anno_sets)):
					list_anno = list(anno_sets[j])
					list_anno.sort()
					clust_to_print = [float(c) for c in clusters[int(old_clust)].split(",")]
					scaled_clust_to_print = [scale * x for x in clust_to_print]
					file.write(file_anno + "_" + old_clust + "\t" + ",".join(list_anno) + "\t" + ",".join(map(str, scaled_clust_to_print)) + "\n")
			#Clear and update the annotation list.
			anno_sets.clear()
			anno_sets.append(set(split_anno))
		
		#Update the current cluster.
		old_clust = new_clust
			
	#Close the file.
	file.close()


if __name__ == "__main__":
	main()