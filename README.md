<h1>Getting Started</h1>
    <p>This repository contains the shapes and in-house scripts used in our paper <b>"Regulatory Element Annotation of the Genome from Chromatin Accessibility Signal Shape us-ing Modified Self-Organizing Maps"</b>. The code in this repository can be used to annotate new chromatin accessibility samples, learn new shapes, or append to a set of existing shapes. We have also included code to replicate our results. Our code is designed for use in a Unix environment and can be run using a command-line interface. 
    <h2>Dependencies</h2>
    <ul>
        <li>Python 3. Our scripts are all written to run in Python 3. If using an older version of Python, they may be incompatible.</li>
        <li>bedtools. This can be installed here: https://github.com/arq5x/bedtools2/releases</li>
        <li>The Unix utilities shuf, cut, and awk. These should be available on most Unix systems.</li>
        <li>The Python modules numpy, scipy, and tqdm.</li>
    </ul>
<h1>Downloading Our Data</h1>
    <p> To download the data used in our paper, run <b>common_scripts/download_all_files.sh</b>. Make sure you run this command from the folder where you wish to download the files. This will download all BAM files, peaks, and ChromHMM mnemonics used in our analyses and concatenate them where needed. This script takes no parameters and will name the BAM files according to cell type.</p>
    <h3>Additional Dependencies</h3>
    <ul>
        <li>bamtools. This can be installed here: https://bioconda.github.io/recipes/bamtools/README.html</li>
        <li>wget. This can be installed here: https://ftp.gnu.org/gnu/wget/</li>
    </ul>
<h1>Learning Shapes (VNSSOM)</h1>
    <p>The code you will need for this task is in the folder <b>shape_learning_scripts</b>. Given the input BAM and ChromHMM mnemonic files, follow these steps to learn shapes and associate them with ChromHMM mnemonics.</p>
    <h3>Additional Dependencies</h3>
    <ul>
        <li>bamtools. This can be installed here: https://bioconda.github.io/recipes/bamtools/README.html</li>
        <li>bigWigToWig. This can be installed here: https://www.encodeproject.org/software/bigwigtowig/</li>
        <li>Tensorflow. This can be installed here: https://www.tensorflow.org/install/.</li>
        <li>bamCoverage. This can be installed here: https://deeptools.readthedocs.io/en/develop/content/installation.html</li>
        <li>The python library pysam</li>
        <li>BedOps. This can be installed here: https://bedops.readthedocs.io/en/latest/index.html</li>
    </ul>
    <h3>Steps</h3>
    <ol>
        <li>Download the Kundaje Lab blacklist file from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/</li>
        <li>Run <b>bedtools merge -i hg38.blacklist.bed -o hg38.blacklist_merged.bed</b> to create a merged blacklist.</li>
        <li>Run <b>convert_bam_to_wig.sh</b> to convert the BAM files to chromosome-specific WIG files. This requires the following parameters to be specified:
        <ul>
            <li> <b>-n:</b> The name of the cell line (e.g. Brain). This will be used to locate the input BAM file and to name the output WIG files.</li>
            <li> <b>-d:</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
            <li> <b>-b:</b> The BAM file used as input.</li>
            <li> <b>-i:</b> The bin size used to generate the WIG file (default: 50 bp).</li>
            <li> <b>-s:</b> The file containing a list of chromosome sizes. This is needed for splitting the BAM file by chromosome.</li>
            <li> <b>-l:</b> The blacklist regions to exclude. This is the merged blacklist file you created.</li>
        </ul>
        </li>
        <li> Run <b>create_regions.sh</b> to create training regions. To create training regions from a permuted WIG file, run <b>create_training_regions_permuted.sh</b>. The requires the following parameters to be specified:
        <ul>
            <li><b>-n</b> The name of the cell line (e.g. Brain)</li>
            <li><b>-d</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
            <li><b>-i</b> The bin size used to generate the WIG file (default: 50 bp)</li>
            <li><b>-r</b> The region size used for splitting (default: 4 kbp)</li>
            <li><b>-w</b> The directory containing the WIG file</li>
            <li><b>-o</b> The directory to contain the split regions</li>
            <li><b>-m</b> The margin to use in splitting the regions</li>
            <li><b>-f</b> The factor to use in splitting the regions</li>
        </ul>
        </li>
        <li> Run <b>learn_shapes_for_chrom_vnssom.sh</b> to learn shapes on one chromosome using our method. You may also use the script <b>learn_all_shapes_vnssom_pbs</b> if you wish to learn shapes on all chromosomes and are running on a supercomputer that uses a PBS job system. To learn shapes using other methods, similar scripts are available, with the extensions <b>_cagt</b>, <b>_chromhmmperm</b>, <b>_signal</b>, and <b>_som</b>. This requires the following parameters to be specified:
        <ul>
            <li><b>-d:</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
            <li><b>-h:</b> The ChromHMM file used for intersecting.</li>
            <li><b>-i:</b> The bin size used to generate the WIG file (default: 50 bp)</li>
            <li><b>-r:</b> The size of the input regions (default: 4000)</li>
            <li><b>-t:</b> The cutoff to use for cross-correlation significance.</li>
            <li><b>-a:</b> Directory containing training regions</li>
            <li><b>-u:</b> Percentile cutoff file</li>
            <li><b>-c:</b> The chromosome name</li>
        </ul>
        </li>
    </ol>
    <h3>Description of Helper Scripts</h3>
    <ul>
        <li><b>create_index_pysam.py:</b> Indexes a BAM file using Pysam. This is necessary because, if you index using bamtools, the <b>bamCoverage</b> utility cannot process the file.</li>
        <li><b>split_regions.py:</b> Splits the WIG file into regions to use for training the VNSSOM.</li>
        <li><b>vnssom.py:</b> The central VNSSOM script. It learnes the shapes given the training regions.</li>
        <li><b>som.py:</b> The vanilla version of the central VNSSOM script.</li>
        <li><b>permute_chromhmm.py:</b> Permutes the ChromHMM mnemonics.</li>
        <li><b>permute_wig.py:</b> Permutes the WIG signal intensities.</li>
        <li><b>extract_signal.py:</b> Extracts pickled input regions and stores them in CSV files for use by CAGT.</li>
        <li><b>run_cagt.m:</b> Runs CAGT. Note that MATLAB is required to run this script, because CAGT is implemented in MATLAB. CAGT must also be installed; you can download it here: https://github.com/sofiakp/cagt/tree/master/matlab</li>
        <li><b>merge_shifted.py:</b> Consolidates shifted shapes learned by the VNSSOM using cross-correlation.</li>
        <li><b>make_shape_bed.py:</b> Annotate each region with its closest shape learned by the VNSSOM.</li>
        <li><b>find_chromhmm_distrib.py:</b> Finds the distribution of ChromHMM mnemonics across each shape.</li>
        <li><b>signal_chromhmm_distrib.py:</b> Finds the distribution of ChromHMM mnemonics across each signal intensity.</li>
    </ul>
    <h3>Output</h3>
    <ul>
        <li>A WIG file for each chromosome in the directory <b>wig</b> in the base directory.</li>
        <li>Percentile cutoffs for computing crossing count in the <b>percentile_cutoffs</b> folder of the directory specified.</li>
        <li>The training regions in the <b>shifted</b> directory of the directory specified.</li>
        <li>The shapes learned by the VNSSOM and shifting procedure in the directory <b>vnssom_output_shifted</b> in the base directory.</li>
        <li>The training regions annotated with shapes in the directory <b>vnssom_anno_beds</b> in the base directory.</li>
        <li>The sorted list of regions annotated with these shapes in the directory <b>anno_beds_sorted</b> in the base directory.</li>
        <li>The intersections between our shapes and the ChromHMM regulatory annotations in <b>vnssom_intersects</b> in the base directory.</li>
        <li>The shapes with their associated ChromHMM mnemonics in the directory <b>vnssom_chromhmm_distrib</b> in the base directory.</li>
    </ul>
    

<h1>Replicating Our Results</h1>
<p>To download the data used in our analysis, run the script <b>download_all_files.sh</b> under <b>common_scripts</b>. You will need to run it from the location where you wish to save the BAM files.</p>
<h3>Additional Dependencies</h3>
<p>To download our data, you will need a system with wget (Most Unix systems should have this). Otherwise, you can download the data manually. You will also need the Python packages glob, pandas, sklearn, matplotlib, and seaborn to run the remaining scripts. Please note that all images will be saved to a file; you do not need a graphical user interface to run this code.</p>
<h3>Replicating Baseline Results</h3>
<ul><li>To replicate the permuted ChromHMM experiment, run <b>associate_from_perm_chromhmm.sh</b>. Options are <b>-s</b>, the source cell line from which the shapes were learned, <b>-d</b>, the name of the cell line to annotate, <b>-b</b>, the directory where to save output, <b>-w</b>, the location of the WIG chromosomes from the source cell line, <b>-a</b>, the location of the WIG chromosomes to annotate, <b>-c</b>, the ChromHMM file to permute, and <b>-p</b>, the name of the permuted ChromHMM file to save.</li>
<li>To replicate the permuted WIG signal experiment, use <b>learn_from_perm_wig.sh</b>. The options are <b>-n</b>, the source cell line from which the shapes were learned, <b>-d</b>, the base directory from where files are saved, <b>-c</b>, the file with the ChromHMM annotations, <b>-i</b>, the bin size (default 50 bp), <b>-r</b>, the region size (default 4000 bp), and <b>-w</b>, the WIG file.</li> 
<li>To annotate promoters using transcription start sites, use <b>get_tss_predictions.sh</b>.</li> 
<li>To annotate promoters using RPKM, use <b>predict_from_rpkm.py</b> in the <b>annotation_scripts</b> directory.</li> 
<li>To replicate the PEAS experiment, use <b>associate_non_promoters_peas.sh</b>. Options are the same as for the promoted signal.</li> 
<li>To replicate the CAGT experiment, use <b>cagt_prep.py</b>, the CAGT code from https://code.google.com/archive/p/cagt/source/default/source, and <b>cagt_associate.sh</b>. Note that you will need to have Matlab installed on your system.</li></ul> 
<h3>Replicating Figures</h3>
<ul><li><b>Fig. 3: </b> Run <b>plot_precision_recall.py</b> followed by <b>plot_precision_recall_all.py</b> (for part A) and <b>plot_true_distribs_all.py</b> (for part B).</li>
<li><b>Fig. 4: </b> Run <b>run_all_chromosome_iterations.sh</b> on each cell type to generate the density plot followed by <b>save_precision_recall.py</b> and <b>plot_precision_recall_densities.py</b>.</li>
<li><b>Fig. 5: </b> Run <b>plot_precision_recall_nobaselines.py</b>.</li>
<li><b>Fig. 6: </b> Run <b>plot_precision_recall_nopromoter_abovethreshonly.py</b>.</li>
<li><b>Supplementary Fig. 1: </b> Run <b>plot_wig_distribs_violin.py</b>.</li>
<li><b>Supplementary Fig. 4: </b> Run <b>plot_chromhmm_distribs_violin.py</b>.</li>
<li><b>Supplementary Fig. 5: </b> Run <b>annotation_similarity_heatmap.py</b>.</li>
<li><b>Supplementary Fig. 6: </b> Run <b>print_annotated_shapes.py</b>.</li>
<li><b>Supplementary Fig. 7, 8, and 9: </b> Run <b>plot_precision_recall.py</b>.</li>
<li><b>Supplementary Fig. 10: </b> Run <b>plot_crosscorr_distrib.py</b>.</li>
<li><b>Supplementary Fig. 11: </b> Run <b>plot_precision_recall_nobaselines.py</b>.</li></ul>
<h1>Annotating New Samples</h1>
    <p>The code you will need for this task is in the folder <b>annotation_scripts</b>. You should only need to run <b>do_annotation.sh</b> or, to run the TSS AND SOM-VN model, <b>do_annotation_and.sh</b>. By default, the script is set to run on each of the 24 human chromosomes. You can modify this as needed by changing the value of the <code>CHROMS</code> variable.</p>
    <h3>Options</h3>
        <p><b>-n:</b> The name you wish to assign to your sample. <em>Required</em></p>
        <p><b>-d:</b> The directory where you would like all files and annotations to be saved. <em>Required</em></p>
        <p><b>-s:</b> The path to the shape file to use for annotation. <em>Required</em></p>
        <p><b>-b:</b> The path to the input BAM file. <em>Required</em></p>
        <p><b>-i:</b> The size of bins you wish to use in generating your WIG file (in bp). The default is 50. <em>This must be the same bin size used to learn the shapes.</em></p>
        <p><b>-w:</b> The path to wig_split. <em>Required</em></p>
        <p><b>-r:</b> The size of region you wish to annotate (in bp). The default is 8000. </p>
        <p></p>
        <p><b>Note:</b> If running the TSS AND SOM-VN model, you must first run get_tss_promoters.sh. This script contains the following options:</p>
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