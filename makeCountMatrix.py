#Parallel processing
from joblib import Parallel, delayed
import multiprocessing

#Import numpy for calculations
import numpy as np

#Import glob for obtaining the files.
import glob, os

#Import sys for obtaining command line args.
import sys

"""
For each datum in each of the window size files, match it to its
closest cluster in the cluster files. Print the region, its closest
match, and the ambiguity of the match to a BED file.
"""
def main():

	#Window size list and list of chromosomes are variable.
	#Allow user to choose number of iterations, initial learning
	#rate, and initial neighborhood size.
	#20, 50, 100, 250, 
	input = sys.argv[1]
	clusters = sys.argv[2]
	output = sys.argv[3]
	
	#Open all input files.
	in_file = open(input, 'r')
	cluster_file = open(clusters, 'r')
	
	#Make 2-d matrix of counts.
	cluster_count = get_cluster_count(cluster_file)
	annotation_list = ["1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx", "6_TxWk",\
	"7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF/Rpts", "13_Het", "14_TssBiv",\
	"15_EnhBiv", "16_ReprPC", "17_ReprPCWk", "18_Quies"]
	#annotation_list = ["1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "5_Tx", "6_TxWk",\
	#"7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF/Rpts", "13_Het", "14_TssBiv",\
	#"15_EnhBiv", "16_ReprPC", "17_ReprPCWk", "18_Quies"]
	#annotation_list = ["1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", \
	#"8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies"]
	#rangeOfVals = range(1, 51)
	#annotation_list = ["E" + str(i) for i in rangeOfVals]
	count_matrix = np.zeros((len(annotation_list), cluster_count))
	
	#Tally up the counts for each cluster for each annotation.
	next_line = in_file.readline()
	while(next_line):
	
		#Get the annotation and cluster index.
		split_line = next_line.split("\t")
		anno_index = -1
		#if split_line[8] in annotation_list:
		anno_index = annotation_list.index(split_line[8])
		cluster_index = int(split_line[3])
		#if split_line[8] == "13_ReprPC":
		#	anno_index = 16
		
		#Increment the count matrix at the correct position.
		count_matrix[anno_index, cluster_index] += 1

		#Read the next line in the file.			
		next_line = in_file.readline()
		
	#Print the count matrix to a file.
	np.savetxt(output, count_matrix, fmt = "%d", delimiter = ",")

	
"""
Find the number of clusters given the cluster file.
"""
def get_cluster_count(cluster_file):
	
	#Count number of lines in cluster.
	cluster_count = 0
	next_cluster = cluster_file.readline()
	while next_cluster:
		cluster_count += 1
		next_cluster = cluster_file.readline()
		
	return cluster_count
	

if __name__ == "__main__":
	main()