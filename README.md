# Getting Started
This repository contains the shapes and in-house scripts used in our paper *Self-Organizing Maps with Variable Neighborhoods Facilitate Learning of Chromatin Accessibility Signal Shapes Associated with Regulatory Elements*. The code in this repository can be used to annotate new chromatin accessibility samples or learn new shapes. We have also included code to replicate our results. Our code is designed for use in a Unix environment and can be run using a command-line interface.

## Dependencies
* Python 2.7 or higher. Note that your version of Python must be compiled to use UCS2. You can check this by running the following commands on the Python command line:
 `import sys`
 `print sys.maxunicode`
 If it uses UCS2, it should print *65535*.
* wig_split.py. This can be obtained via the [taolib](https://github.com/taoliu/taolib) repository. **Note:** We plan to release a future version without this dependency.
* The GCC compiler. Some our file I/O scripts are written in C and will need to be compiled using GCC. This should be available on most Unix systems.
* bedtools. This can be installed [here](https://github.com/arq5x/bedtools2/releases).
* The Unix utilities shuf, cut, and awk. These should be available on most Unix systems.
* The Python modules numpy, scipy, HTSeq, argparse, logging, sys, itertools, and pysam.
* [Python code to compute the gap statistic](https://github.com/minddrummer/gap/blob/master/gap/gap.py)

* If you plan to learn new shapes, the Python module tensorflow.

## Installing gosr
While a [pre-existing version of gosr exists](https://github.com/wresch/gosr), our framework uses a modified version that preserves zeros in the signal profile. gosr should automatically be downloaded when you download our repository. To install it, simply untar the file file `gosr.tar` provided in this repository.  

## Adding programs to PATH and PYTHONPATH
You will need to add the following paths to your `PYTHONPATH` variable:
* *<installation_path>/SOM-VN/gosr/gosr/tools*
* *<path_to_gap_statistic_code>*

You will need to add the following paths to your PATH variable:
* *<installation_path>/SOM-VN/gosr/bin*
* The directory where you installed bedtools.
* The directory containing the *taolib* folder for `wig_split.py`.

# Obtaining regions to annotate from a BAM file
To segment regions for learning shapes or annotating shapes, you will need to run the following:
1. `gosr binbam -f 0 -n 1000 -t <name_of_sample> <bam_file> 50 <name_of_sample> > <name_of_WIG_file>`
1. `sed '3d' <name_of_WIG_file> > <name_of_WIG_file>_noheader.wig`
1. `mkdir <chrom_wig_dir>`
1. `python wig_split.py <name_of_WIG_file>_noheader.wig <chrom_wig_dir>/<name_of_sample>`
1. `gcc -pthread -lm -o run_get_data ../common_scripts/get_file_data.c`

# Annotating New Samples
You may annotate new samples using the shape files we provide or using new shapes that you have learned. In the pipeline below, *file_containing_shapes* refers to the shape file you are using.
The code you will need for this task is in the folder *annotation_scripts*. To annotate the new samples, you will need to do the following for each chromosome. Here, the *annotation_file* is the location where you wish to store the annotated regions.

1. `python make_annotated_bed.py <regions_to_annotate> <file_containing_shapes> <annotation_file> <chrom_wig_file>  0.0`
 This will generate two files: one that contains the annotated regions and their scores, and one that also includes the annotated regions with the original signal, saved as *<annotation_file>clust*.
1. `bedtools sort -i  <annotation_file> > <annotation_sorted_file>`
1. `bedtools sort -i  <annotation_file>clust > <annotation_sorted_file>clust`
1. `python common_scripts/consolidate_bed.py <annotation_sorted_file> <annotation_consolidated_file>`

# Learning Shapes
The code you will need for this task is in the folder *shape_learning_scripts*. To learn shapes for one chromosome, you will need to do the following:
1. `shuf <regions_to_train> > <shuffled_regions_to_train>`
1. `python shift_input.py <shuffled_regions_to_train> <shuffled_regions_to_train> 50 4000 <chrom_wig_file>.chr<chrom> false 0`
1. `python som_vn.py <shuffled_regions_to_train> <som_shapes_learned> <chrom_wig_file>.chr<chrom> 4000 50 0 False`
1. `python remove_by_cutoff.py <som_shapes_learned> 1 <som_shapes_filtered>`
1. `python merge_shifted.py <som_shapes_filtered> <som_shapes_shifted> 0`
1. `python kmeans_shapes.py <som_shapes_shifted> <som_shapes_final>`
1. `python make_shape_bed.py <regions_to_annotate> <som_shapes_final> <regions_annotated> 0
1. `bedtools sort -i  <regions_annotated> > <regions_annotated_sorted>`
1. `python consolidate.py <regions_annotated_sorted> <regions_annotated_final>`
1. `cut -d$'\t' -f 1,2,3,4,5 <regions_annotated_final> > <regions_annotated_final>.bed`
1. `awk '{ print $6}' <regions_annotated_final> > <clusters_annotated_final>`
1. `cut -d$'\t' -f 7,8,9,10 <regions_annotated_final> > <scores_annotated_final>`
1. `bedtools intersect -wao -a <regions_annotated_final>.bed -b <chromhmm_mnemonic_file> > <intersects>`
1. `bedtools sort -i <intersects> > <intersects_sorted>`
1. `python consolidate_chromHMM.py <intersects_sorted> <som_shapes_final> <shapes> <chrom_wig_file>.chr<chrom> <name_of_sample> <regions_to_annotate> 0`

To merge shapes learned across multiple chromosomes, run the following:
`python common_scripts/merge_significant.py <shapes> <shapes_all> <shapes_log>`

# Replicating Our Results

## Downloading the data
All of our data can be downloaded by running `download_all_files.sh` in the location where you wish to save your BAM files. They will be index by the cell type and SPOT score.

## Evaluating null models
To evaluate null models, you will need to learn shapes using SOM-VN with unpermuted input and SOM-VN with permuted input. You will then need to associate the learned shapes with both orignal, unpermuted ChromHMM RE and permuted ChromHMM RE.

### Learning SOM-VN shapes
1. Run the following scripts. Note that they are designed to be run in a SLURM supercomputing environment, but can be modified to be run on a local machine.
  * `learn_shapes_A549_slurm.sh`
  * `learn_shapes_b_cell_low_slurm.sh`
  * `learn_shapes_b_cell_high_slurm.sh`
  * `learn_shapes_Brain_low_slurm.sh`
  * `learn_shapes_Brain_high_slurm.sh`
  * `learn_shapes_GM12878_slurm.sh`
  * `learn_shapes_H1_low_slurm.sh`
  * `learn_shapes_H1_high_slurm.sh`
  * `learn_shapes_HeLa_low_slurm.sh`
  * `learn_shapes_HeLa_high_slurm.sh`
  * `learn_shapes_heart_low_slurm.sh`
  * `learn_shapes_heart_high_slurm.sh`
  * `learn_shapes_stomach_low_slurm.sh`
  * `learn_shapes_stomach_high_slurm.sh`
2. Merge shapes learned across chromosomes. To do this, run:
	`python common_scripts/merge_significant.py $BASE_PATH/anno_A549_consolidated <shapes_all_A549> <shapes_log_A549>`
	...
	`python common_scripts/merge_significant.py $BASE_PATH/anno_stomach_high_consolidated <shapes_all_stomach_high> <shapes_log_stomach_high>`
	
### Learning SOM-VN shapes from permuted signal
1. Run the following scripts. Note that they are designed to be run in a SLURM supercomputing environment, but can be modified to be run on a local machine.
  * `learn_from_perm_wig_A549_slurm.sh`
  * `learn_from_perm_wig_b_cell_low_slurm.sh`
  * `learn_from_perm_wig_b_cell_high_slurm.sh`
  * `learn_from_perm_wig_Brain_low_slurm.sh`
  * `learn_from_perm_wig_Brain_high_slurm.sh`
  * `learn_from_perm_wig_GM12878_slurm.sh`
  * `learn_from_perm_wig_H1_low_slurm.sh`
  * `learn_from_perm_wig_H1_high_slurm.sh`
  * `learn_from_perm_wig_HeLa_low_slurm.sh`
  * `learn_from_perm_wig_HeLa_high_slurm.sh`
  * `learn_from_perm_wig_heart_low_slurm.sh`
  * `learn_from_perm_wig_heart_high_slurm.sh`
  * `learn_from_perm_wig_stomach_low_slurm.sh`
  * `learn_from_perm_wig_stomach_high_slurm.sh`
2. Merge shapes learned across chromosomes. To do this, run:
	`python common_scripts/merge_significant.py $BASE_PATH/anno_A549_consolidated_perm <shapes_all_A549_perm> <shapes_log_A549_perm>`
	...
	`python common_scripts/merge_significant.py $BASE_PATH/anno_stomach_high_consolidated_perm <shapes_all_stomach_high_perm> <shapes_log_stomach_high_perm>`
	
### Associating SOM-VN shapes with permuted RE
1. Run `python permute_chromhmm.py <chromhmm_mnemonic_file> <permuted chromhmm_mnemonic_file>`
1. Run the following scripts. Note that they are designed to be run in a SLURM supercomputing environment, but can be modified to be run on a local machine.
  * `associate_from_perm_chromhmm_A549_slurm.sh`
  * `associate_from_perm_chromhmm_b_cell_low_slurm.sh`
  * `associate_from_perm_chromhmm_b_cell_high_slurm.sh`
  * `associate_from_perm_chromhmm_Brain_low_slurm.sh`
  * `associate_from_perm_chromhmm_Brain_high_slurm.sh`
  * `associate_from_perm_chromhmm_GM12878_slurm.sh`
  * `associate_from_perm_chromhmm_H1_low_slurm.sh`
  * `associate_from_perm_chromhmm_H1_high_slurm.sh`
  * `associate_from_perm_chromhmm_HeLa_low_slurm.sh`
  * `associate_from_perm_chromhmm_HeLa_high_slurm.sh`
  * `associate_from_perm_chromhmm_heart_low_slurm.sh`
  * `associate_from_perm_chromhmm_heart_high_slurm.sh`
  * `associate_from_perm_chromhmm_stomach_low_slurm.sh`
  * `associate_from_perm_chromhmm_stomach_high_slurm.sh`
2. Merge shapes learned across chromosomes. To do this, run:
	`python common_scripts/merge_significant.py $BASE_PATH/anno_A549_consolidated_chromhmm_perm <shapes_all_A549_chromhmm_perm> <shapes_log_A549_chromhmm_perm>`
	...
	`python common_scripts/merge_significant.py $BASE_PATH/anno_stomach_high_consolidated_chromhmm_perm <shapes_all_stomach_high_chromhmm_perm> <shapes_log_stomach_high_chromhmm_perm>`
	
### Plotting cross-correlation
1. Run the following scripts:
  * `run_crosscorr_A549_slurm.sh`
  * `run_crosscorr_b_cell_low_slurm.sh`
  * `run_crosscorr_b_cell_high_slurm.sh`
  * `run_crosscorr_Brain_low_slurm.sh`
  * `run_crosscorr_Brain_high_slurm.sh`
  * `run_crosscorr_GM12878_slurm.sh`
  * `run_crosscorr_H1_low_slurm.sh`
  * `run_crosscorr_H1_high_slurm.sh`
  * `run_crosscorr_HeLa_low_slurm.sh`
  * `run_crosscorr_HeLa_high_slurm.sh`
  * `run_crosscorr_heart_low_slurm.sh`
  * `run_crosscorr_heart_high_slurm.sh`
  * `run_crosscorr_stomach_low_slurm.sh`
  * `run_crosscorr_stomach_high_slurm.sh`
1. Make the plots. To do this, run:
	`python meta_analysis_scripts/plot_crosscorr_distrib.py $BASE_PATH/crosscorr_A549/anno $BASE_PATH/crosscorr_A549_perm/anno <name_of_plot_A549> A549`
	...
	`python meta_analysis_scripts/plot_crosscorr_distrib.py $BASE_PATH/crosscorr_stomach_high/anno $BASE_PATH/crosscorr_stomach_high_perm/anno <name_of_plot_stomach_high> "Stomach (High)"`
1. Run hypothesis testing, i.e. `python crosscorr_hypothesis_tests.py $BASE_PATH/crosscorr_A549/anno $BASE_PATH/crosscorr_A549_perm/anno`.

## Evaluating SOM-VN against CAGT
To run CAGT, note that you must have MATLAB installed. In addition, while [CAGT is available for download separately](https://github.com/sofiakp/cagt/tree/master/matlab), our local copy of CAGT contains a slight modification that prevents the program from stopping early in the case of slow convergence.
1. Learn SOM-VN shapes as in (###Learning SOM-VN shapes from permuted signal).
1. Run the following commands before starting CAGT.
  1. `mkdir /data/eichertd/som_vn_data/matlab_matrix_<cell_type>`
  1. `mkdir /data/eichertd/som_vn_data/cagt_out_GM12878/`
1. For each chromosome and cell type, run the following:
  1. `matlab -nodisplay -nodesktop -r "run_cagt('/data/eichertd/som_vn_data/training_<cell_type>_shifted/chrom<chrom>','/data/eichertd/som_vn_data/matlab_matrix_<cell_type>','<chrom>','/data/eichertd/som_vn_data/cagt_out_<cell_type>/<chrom>.csv','$BASE_Path/cagt/trunk/matlab/src/', '0.8', '1000', '40')"`
1. Run the following scripts:
  * `run_cagt_analysis_A549_slurm.sh`
  * `run_cagt_analysis_b_cell_low_slurm.sh`
  * `run_cagt_analysis_b_cell_high_slurm.sh`
  * `run_cagt_analysis_Brain_low_slurm.sh`
  * `run_cagt_analysis_Brain_high_slurm.sh`
  * `run_cagt_analysis_GM12878_slurm.sh`
  * `run_cagt_analysis_H1_low_slurm.sh`
  * `run_cagt_analysis_H1_high_slurm.sh`
  * `run_cagt_analysis_HeLa_low_slurm.sh`
  * `run_cagt_analysis_HeLa_high_slurm.sh`
  * `run_cagt_analysis_heart_low_slurm.sh`
  * `run_cagt_analysis_heart_high_slurm.sh`
  * `run_cagt_analysis_stomach_low_slurm.sh`
  * `run_cagt_analysis_stomach_high_slurm.sh`
1. Merge shapes learned using CAGT for each cell type, e.g.
  `python ../common_scripts/merge_significant.py $BASE_PATH/anno_A549_consolidated_cagt <shapes_all_A549_cagt> <shapes_log_A549_cagt>`
1. Run `python plot_precision_recall_nobaselines.py $BASE_PATH/annotated_merged_<cell_type>/ $BASE_PATH/annotated_consolidated_<cell_type>/ $BASE_PATH/pr_<cell_type>.png "<cell_type> SOM-VN" $BASE_PATH/wig/<cell_type>/<cell_type>.chr` to plot precision and recall for SOM-VN.
1. Run `python meta_analysis_scripts/save_precision_recall.py $BASE_PATH/annotated_merged_<cell_type>/ $BASE_PATH/annotated_consolidated_<cell_type>/ <wig_chrom_dir> <cell_type_precision_recall_file>`
1. Run `python plot_precision_recall_nobaselines.py $BASE_PATH/annotated_merged_<cell_type>_cagt/ $BASE_PATH/annotated_consolidated_<cell_type>_cagt/ $BASE_PATH/pr_<cell_type>_cagt.png "<cell_type> CAGT" $BASE_PATH/wig/<cell_type>/<cell_type>.chr` to plot precision and recall for CAGT.
1. Run `python meta_analysis_scripts/save_precision_recall.py $BASE_PATH/annotated_merged_<cell_type>_cagt/ $BASE_PATH/annotated_consolidated_<cell_type>_cagt/ <wig_chrom_dir> <cell_type_precision_recall_file>`

## Evaluating SOM-VN above a threshold
Run `python meta_analysis_scripts/save_precision_recall_threshold.py $BASE_PATH/annotated_merged_<cell_type>/ $BASE_PATH/annotated_consolidated_<cell_type>/ <wig_chrom_dir> <cell_type_precision_recall_file>`

## Evaluating cross-cell-type performance
1. Run the following scripts:
  * `annotate_b_cell_low_from_high.sh`
  * `annotate_HeLa_low_from_high.sh`
  * `annotate_Heart_low_from_high.sh`
  * `annotate_Stomach_low_from_high.sh`
1. Run `mkdir <b_cell_precision_recall_directory>`.
1. Run `python meta_analysis_scripts/save_precision_recall.py $BASE_PATH/annotated_merged_b_cell_high_b_cell_low/ $BASE_PATH/annotated_consolidated_b_cell_high_b_cell_low/ <b_cell_wig_chrom_dir> <b_cell_precision_recall_directory>` 

## Evaluating cross-chromosome performance

## Evaluating performance using magnitude only.
1. Run `shape_learning_scripts/magnitude_all.sh` to associate magnitude with RE on each chromosome of each cell type.
1. Run `python merge_significant_magnitude.py $BASE_PATH/anno_<cell_type>_magnitude_consolidated/ $BASE_PATH/magnitudes_<cell_type> $BASE_PATH/magnitudes_<cell_type>.log`
1. Run `annotation_scripts/annotate_magnitude_all.sh` to annotate each chromosome of each cell type with the learned magnitude associations for that cell type.
1. Run `python meta_analysis_scripts/save_precision_recall.py $BASE_PATH/annotated_merged_magnitude_<cell_type>/ $BASE_PATH/annotated_consolidated_magnitude_<cell_type>/ <wig_chrom_dir> <cell_type_precision_recall_file>`

## Computing coverage
For each sample, run the following:
1. For each chromosome and RE, run `awk '($4 == "<RE>")' $BASE_PATH/annotated_consolidated_<sample>/anno<chrom>.bed > $BASE_PATH/annotated_consolidated_<sample>/anno<chrom>_known_<RE>.bed`. A shortcut is to run `./get_re_specific_anno.sh`.
1. `cat $BASE_PATH/annotated_consolidated_<sample>/anno*_known_<RE>.bed > $BASE_PATH/annotated_consolidated_<sample>/anno_all_known_<RE>.bed`, and similarly for additional samples.
1. `bedtools intersect -a $BASE_PATH/annotated_consolidated_<sample>/anno_all_known_<RE>.bed -b $PEAKS/<sample>.bed -u > annotations_<RE>_<cell_type>_peak.bed`
1. `wc -l $BASE_PATH/annotated_consolidated_<sample>/anno_all_known_<RE>.bed`
1. `wc -l annotations_<RE>_<cell_type>_peak.bed`

## Computing ChromHMM percentages
Run the following for each cell type.
1. Run `cat $CHROMHMM/b_cell_E032.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'` to compute the total coverage of the entire ChromHMM file.
1. Run `awk '($4 == "6_EnhG" || $4 == "7_Enh" || $4 == "12_EnhBiv")' $CHROMHMM/b_cell_E032.bed > $CHROMHMM/b_cell_E032_enhancer.bed` to obtain a ChromHMM file with only enhancers.
1. Run `awk '($4 == "1_TssA" || $4 == "2_TssAFlnk" || $4 == "10_TssBiv" || $4 == "11_BivFlnk")' $CHROMHMM/b_cell_E032.bed > $CHROMHMM/b_cell_E032_promoter.bed` to obtain a ChromHMM file with only promoters.
1. Run `awk '($4 == "9_Het" || $4 == "15_Quies")' $CHROMHMM/b_cell_E032.bed > $CHROMHMM/b_cell_E032_weak.bed` to obtain a ChromHMM file with only weak RE.
1. Run `cat $CHROMHMM/b_cell_E032_<RE>.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'` to compute the total coverage of ChromHMM for each RE (promoter, enhancer, weak).

## Compute overlap between results for single replicates.
First, download and unzip the peaks for each high quality brain tissue replicate.
1. `wget https://www.encodeproject.org/files/ENCFF018ATG/@@download/ENCFF018ATG.bed.gz`
1. `gunzip ENCFF018ATG.bed.gz`
1. `mv ENCFF018ATG.bed Brain_rep1_peaks.bed`
1. `wget https://www.encodeproject.org/files/ENCFF936VAD/@@download/ENCFF936VAD.bed.gz`
1. `gunzip ENCFF936VAD.bed.gz`
1. `mv ENCFF936VAD.bed Brain_rep2_peaks.bed`
Next, run IDR:
1. `idr --samples Brain_rep{1,2}_peaks.bed`
Follow the same steps for A549 and GM12878, where the peak files to download are as follows:
* A549 Rep 1: *https://www.encodeproject.org/files/ENCFF529HMB/@@download/ENCFF529HMB.bed.gz*
* A549 Rep 2: *https://www.encodeproject.org/files/ENCFF823UOG/@@download/ENCFF823UOG.bed.gz*
* GM12878 Rep 1: *https://www.encodeproject.org/files/ENCFF598KWZ/@@download/ENCFF598KWZ.bed.gz*
* GM12878 Rep 2: *https://www.encodeproject.org/files/ENCFF073ORT/@@download/ENCFF073ORT.bed.gz*
Run the following files to learn shapes for each replicate:
1. `learn_shapes_Brain_rep1_slurm.sh`
1. `learn_shapes_Brain_rep2_slurm.sh`
1. `learn_shapes_GM12878_rep1_slurm.sh`
1. `learn_shapes_GM12878_rep2_slurm.sh`
1. `learn_shapes_A549_rep1_slurm.sh`
1. `learn_shapes_A549_rep2_slurm.sh`
1. For each cell type, `python common_scripts/merge_significant.py $BASE_PATH/anno_<sample_name>_consolidated_perm <shapes> <shapes_log>`
Run the following files to annotate each cell type for each replicate:
1. `annotate_Brain_rep1.sh`
1. `annotate_Brain_rep2.sh`
1. `annotate_GM12878_rep1.sh`
1. `annotate_GM12878_rep2.sh`
1. `annotate_A549_rep1.sh`
1. `annotate_A549_rep2.sh`
Prepare chromatin state assignments as follows for each sample:
1. `cd $BASE_PATH/annotated_consolidated_<sample_name>/`
1. `rm anno*clust.bed`
1. `cp anno*.bed > $BASE_PATH/anno_all.bed
1. sort -k1,1 -k2,2n $BASE_PATH/annotated_consolidated_<sample_name>/all_anno.bed > $BASE_PATH/annotated_consolidated_<sample_name>/all_anno_sorted.bed
1. `awk '($4 == "Enhancer")' $BASE_PATH/all_anno_sorted.bed > $BASE_PATH/all_anno_sorted_enhancer.bed`
1. `awk '($4 == "Promoter")' $BASE_PATH/all_anno_sorted.bed > $BASE_PATH/all_anno_sorted_promoter.bed`
1. `awk '($4 == "Weak")' $BASE_PATH/all_anno_sorted.bed > $BASE_PATH/all_anno_sorted_weak.bed`
For each cell type and RE, compute the overlap between replicates:
1. `bedtools intersect -a annotated_consolidated_<cell_type>_rep1/all_anno_sorted_<RE>.bed -b annotated_consolidated_<cell_type>_rep2/all_anno_sorted_<RE>.bed -u > annotations_<RE>_<cell_type>_a.bed`
1. `bedtools intersect -a annotated_consolidated_<cell_type>_rep2/all_anno_sorted_<RE>.bed -b annotated_consolidated_<cell_type>_rep1/all_anno_sorted_<RE>.bed -u > annotations_<RE>_<cell_type>_b.bed`
1. `wc -l annotated_consolidated_<cell_type>_rep1/all_anno_sorted_<RE>.bed`
1. `wc -l annotations_<RE>_<cell_type>_a.bed`
1. `wc -l annotated_consolidated_<cell_type>_rep2/all_anno_sorted_<RE>.bed`
1. `wc -l annotations_<RE>_<cell_type>_b.bed`



