<h1>Getting Started</h1>
<p>This repository contains the shapess and in-house scripts used in our paper <b>"SOM-VN: Detection and Annotation of Regulatory Elements from Chromatin Accessibility Signal Shape"</b>. The code in this repository can be used to annotate new chromatin accessibility samples, learn new shapess or append to our existing shapess, and replicate our results. It is designed for use in a Unix environment and can all be run using a command-line interface. This README has been divided into each of these three subtasks for your convenience. 
 <h2>Dependencies</h2>
<ul><li>Python 2. Our script generates a WIG file using the tool <b>gosr binbam</b>, for which we have included a modified version. This tool is only compatible with Python 2.</li>
<li>wig_split.py. This can be obtained via the following repository: https://github.com/taoliu/taolib</li>
<li>The GCC compiler. A couple of our file I/O scripts are written in C and will need to be compiled using GCC. This should be available on most Unix systems.</li>
 <li>bedtools. This can be installed here: https://github.com/arq5x/bedtools2/releases</li>
  <li>The Unix utilities shuf, cut, and awk. These should be available on most Unix systems.</li>
<li>The Python modules numpy and scipy.</li>
<li>The Python modules HTSeq, argparse, logging, sys, itertools, and pysam. These are necessary for running gosr.</li>
<li>Your version of Python must be compiled to use UCS2. You can check this by running the following commands on the Python command line:
 <p><code>>import sys</code></p>
 <p><code>>print sys.maxunicode</code></p>
 If it uses UCS2, it should print "65535".</li></ul>
 <h2>Installing gosr</h2>
 <p>While a pre-existing version of gosr exists at https://github.com/wresch/gosr, our framework uses a modified version that preserves zeros in the signal profile. gosr should automatically be downloaded when you download our repository. To install it, simply untar the file file <b>gosr.tar</b> provided in this repository.
 <ol></ol>
<h2>Add to System Path</h2>
<ul><li>The directory where you installed bedtools.</li>
<li>The gosr/bin directory under the directory where you installed gosr.</li>
<li>The directory containing the python packages for gosr.</li>
<li>The directory containing wig_split.py.</li></ul>
<h1>Annotating New Samples</h1>
<p>The code you will need for this task is in the folder <b>annotation_scripts</b>. You should only need to run <b>do_annotation.sh</b> or, to run the TSS AND model, <b>do_annotation_and.sh</b>. By default, the script is set to run on each of the 24 human chromosomes. You can modify this as needed by changing the value of the <code>CHROMS</code> variable.</p>
<h3>Options</h3>
<p><b>-n:</b> The name you wish to assign to your sample. <em>Required</em></p>
<p><b>-d:</b> The directory where you would like all files and annotations to be saved. <em>Required</em></p>
<p><b>-s:</b> The path to the shape file to use for annotation. <em>Required</em></p>
<p><b>-b:</b> The path to the input BAM file. <em>Required</em></p>
<p><b>-i:</b> The size of bins you wish to use in generating your WIG file (in bp). The default is 50. <em>This must be the same bin size used to learn the shapes.</em></p>
<p><b>-w:</b> The path to wig_split. <em>Required</em></p>
<p><b>-r:</b> The size of region you wish to annotate (in bp). The default is 8000. </p>
<p></p>
<p><b>Note:</b> If running the TSS AND model, you must first run get_tss_promoters.sh. This script contains the following options:</p>
<p><b>-n:</b> The name you wish to assign to your sample. <em>Required</em></p>
<p><b>-d:</b> The directory where you would like all files and annotations to be saved. <em>Required</em></p>
<p><b>-t:</b> The path to a file containing the transcription start sites.  <em>Required</em></p>
<p><b>-g:</b> The path to a file containing the chromosome sizes for the genome.  <em>Required</em></p>
<h3>Output</h3>
<ul><li>A WIG file for all chromosomes in the BAM file. This is in the base directory and starts with the name of the cell line.</li>
<li>A WIG file for each chromosome in the directory <b>annotation_files</b>.</li>
<li>The signal data for each region to annotate in the directory <b>wig_chroms_annotation</b>.</li>
<li>An annotated file for each chromosome with the signal data for each region in the directory <b>annotated</b>.</li>
<li>A sorted annotated file with the annotations and signal in the directory <b>annotated_sorted</b>.</li>
<li>An annotated file without overlap and including the signal in the directory <b>annotated_consolidated</b>.</li>
<li>A BED file of the annotations in the directory, excluding the signal, in the directory <b>annotated_consolidated</b>.</li>
<li>A file with only the signal data for each annotation in the directory <b>annotated_consolidated</b>.</li></ul>
<h3>Description of Helper Scripts</h3>
<ul><li><b>make_annotated_bed.py:</b> Annotates a chromosome given the shapes, the signals for each region, and the WIG file.</li>
<li><b>consolidate_bed.py:</b> Consolidates a set of sorted annotations to remove overlap using the ambiguity score described in the paper.</li>
<li><b>get_file_data.c:</b> Segments WIG file into regions for annotation.</li></ul>
<h1>Creating / Appending to a Shape List</h1>
<p>The code you will need for this task is in the folder <b>shape_learning_scripts</b>. You should only need to run <b>learn_shapes.sh</b>. By default, the script is set to run on each of the 24 human chromosomes. You can modify this as needed by changing the value of the <code>CHROMS</code> variable.</p>
<h3>Additional Dependencies</h3>
This code requires the Tensorflow framework, which can be installed here: https://www.tensorflow.org/install/. It also requires the imp module to import the gap statistic code from an external repository. Alternatively, you can download the gap statistic code here: https://github.com/minddrummer/gap.
<h3>Options</h3>
<p><b>-n:</b> The name you wish to assign to your training sample. <em>Required</em></p>
<p><b>-d:</b> The directory where you would like your shapes and all intermediary files to be saved. <em>Required</em></p>
<p><b>-b:</b> The path to the ChromHMM regions. <em>Required</em></p>
<p><b>-b:</b> The path to the training BAM file. <em>Required</em></p>
<p><b>-i:</b> The size of bins you wish to use in training (in bp). The default is 50. <em>This must be the same bin size as you wish to use for annotation.</em></p>
<p><b>-w:</b> The path to wig_split. <em>Required</em></p>
<p><b>-r:</b> The size of region you wish to annotate (in bp). The default is 4000. </p>
<p><b>-s:</b> The name of the shape list to save. To append to one of our shape lists, download it and use that name here. <em>Required</em></p>

<h3>Output</h3>
<ul><li>A WIG file for all chromosomes in the BAM file. This is in the base directory and starts with the name of the cell line.</li>
<li>A WIG file for each chromosome in the directory <b>wig_chroms</b>.</li>
<li>The signal data for each region to use in training in the directory <b>training</b>.</li>
 <li>The signal data for each region to annotate with shapes during training in the directory <b>training_anno</b>.</li>
 <li>The signal data for each sub-region found according to the Region Segmentation procedure in the paper in the directory <b>training_shifted</b>.</li>
 <li>The shapes learned by the SOM-VN in the directory <b>som_output</b>.</li>
 <li>The shapes with at least one mapping in the last iteration in the directory <b>som_output_filtered</b>.</li>
 <li>The shapes with all shifted regions merged in the directory <b>som_output_shifted</b>.</li>
 <li>The shapes clustered using K-means in the directory <b>som_output_final</b>.</li>
 <li>The regions in <b>training_anno</b> annotated with shapes from <b>som_output_final</b> in the directory <b>anno_beds</b>.</li>
<li>The sorted list of regions annotated with these shapes in the directory <b>anno_beds_sorted</b>.</li>
<li>The consolidated non-overlapping list of regions annotated with these shapes in the directory <b>anno_beds_final</b>.</li>
<li>The intersections between our shapes and the ChromHMM regulatory annotations in <b>anno_intersects</b>.</li>
<li>The shapes containing all learned shapes on all chromosomes in the file <b>shapes_all</b>.</li>
<li>The shapes containing a consolidated list of shapes merged using cross-correlation in the file you specified.</li>
<li>Figures corresponding to the breakdown of ChromHMM annotations in each of our shapes in the directory <b>annotation_distribution_figs</b>.</li></ul>
<h3>Description of Helper Scripts</h3>
<ul><li><b>get_file_data.c:</b> Generates sub-regions to use both in training and in finding the associations between shapes and ChromHMM annotations.</li>
 <li><b>shift_input.py:</b> Collects sub-regions to use in training the SOM-VN given the training regions with margins.</li>
 <li><b>som_vn.py:</b> The central SOM-VN script. It learnes the shapes given the input regions.</li>
 <li><b>remove_by_cutoff.py:</b> Removes shapes learned by the SOM-VN that did not have any regions mapping to them in the last iteration of the algorithm.</li>
 <li><b>merge_shifted.py:</b> Consolidates shifted shapes learned by the SOM-VN using cross-correlation.</li>
 <li><b>kmeans_shapes.py:</b> Clusters shapes learned by the SOM-VN using K-means.</li>
 <li><b>make_shape_bed.py:</b> Annotate each region with its closest shape learned by the SOM-VN.</li>
 <li><b>consolidate_chromHMM.py:</b> Creates a shapes of regulatory-associated shapes given the intersections between our shapes and ChromHMM annotations.</li></ul>
<h1>Replicating Our Results</h1>
<p>To download the data used in our analysis, run the script <b>download_all_files.sh</b>. You will need to run it from the location where you wish to save the BAM files.</p>
<p>To generate a shapes and annotate with permuted ChromHMM annotations, run <b>associate_from_perm_chromhmm.sh</b>. Options are the <b>-s</b>, the source cell line from which the shapes were learned, <b>-d</b>, the name of the cell line to annotate, <b>-b</b>, the directory where to save output, <b>-w</b>, the location of the WIG chromosomes from the source cell line, <b>-a</b>, the location of the WIG chromosomes to annotate, <b>-c</b>, the ChromHMM file to permute, and <b>-p</b>, the name of the permuted ChromHMM file to save.</p>
<p>To learn chromosome-specific models on permuted WIG data, use <b>learn_from_perm_wig.sh</b>. The options are <b>--wig</b>, the WIG file to permute and <b>--base_dir</b>, the directory where to save output.</p> 
<p>To annotate promoters using transcription start sites, use <b>get_tss_predictions.sh</b>. The options are <b>--tss</b>, the file containing transcription start sites, <b>--genome</b>, the genome used for alignment of the input BAM file, and <b>--base_dir</b>, the directory where to save output.</p> 
<p>To analyze and plot the precision and recall, use <b>do_precision_recall_analysis.sh</b>. By default, this runs on all cell lines tested in our study. To change this functionality, you should modify the bash script.</p>
<p>To generate violin plots, averaged precision and recall plots and averaged distribution plots over all chromosomes and cell types, and shapes with annotation and TSS counts, use <b>make_plots_all_cells.sh</b>. You can modify input directories at the beginning of the script.</p>
<p>To generate the annotation similarity heatmap over multiple window sizes, use <b>make_cell_specific_plots.sh</b>. You can modify input directories at the beginning of the script.</p>
<h3>Additional Dependencies</h3>
<p>To download our data, you will need a system with wget (Most Unix systems should have this). Otherwise, you can download the data manually. You will also need the Python packages glob, pandas, sklearn, matplotlib, and seaborn to run the remaining scripts. Please note that all images will be saved to a file; you do not need a graphical user interface to run this code.</p>
<h3>Description of Helper Scripts</h3>
<ul><li><b>plot_wig_distribs.py:</b> Plots the distribution of RPKM given three genome-wide WIG files (i.e. Figure 3) and saves them to the file <b>wig_distribs.png</b>.</li>
 <li><b>plot_annotation_distribs.py:</b> Plots the distribution in width of ChromHMM annotations (i.e. Figure 7) and saves them to the file <b>all_annotation_distribs.png</b>.</li>
  <li><b>compute_validity.py:</b> Given the set of shapes learned for all window sizes, computes the Davies-Bouldin indices for each chromosome (i.e. Figure 8) and saves them to the heatmap <b>cluster_validity.png</b>.</li>
 <li><b>annotation_similarity_heatmap.py:</b> Given a shapes with all shapes from all chromosomes, plots a heatmap of shape matches by cross-correlation and a distribution of the resulting association types (i.e. Figures 10 and 11) and saves them to the directory <b>annotation_figs</b>.</li>
<li><b>combine_prediction_beds.py:</b> Combines the TSS-based and shape-based predictions and save them to the directory <b>annotated_combined</b>.</li>
<li><b>runGetRegions.c:</b> Given a BED file in sorted order with no overlap, retrieves the associated signals and saves them to a CSV file.</li>
<li><b>plot_precision_recall.py:</b> Computes the precision and recall for all chromosomes in a cell line over all annotation types and over shape-based, TSS-based, TSS + Shape, and permuted categories and saves them to a scatterplot in <b>precision_recall.png</b>.</li><li><b>print_annotated_clusters.py:</b> Plots the annotation count and average TSS count for each cell type annotated by each cluster in a shapes and saves the plot in <b>cluster_plots.png</b>.</li>
<li><b>line_plot_misclassified_distribs.py:</b> Plots the distribution of unknown and polycomb regions and saves the result in the directory <b>misclassified</b>.</li>
<li><b>getBedRegions.c:</b> Outputs signal corresponding to regions in a BED file, filtering out regions above a given window size. This is used in evaluating TSS-based annotation.</li>
<li><b>getBedLines.c:</b> Outputs signal corresponding to regions in a BED file. This is used in evaluating combined TSS + Shape annotations.</li></ul>
