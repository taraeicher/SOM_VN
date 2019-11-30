* [Getting Started](#Getting-Started)
* [Downloading Our Data](#Downloading-Our-Data)
* [Learning Shapes](#Learning-Shapes)
* [Annotating New Samples](#Annotating-New-Samples)
* [Evaluating Ground Truth on Peaks](#Evaluating-Ground-Truth-On-Peaks)
* [Replicating Our Experiments](#Replicating-Our-Experiments)
* [Replicating Our Figures](#Replicating-Our-Figures)

# Getting Started
This repository contains the shapes and in-house scripts used in our paper *Regulatory Element Annotation of the Genome from Chromatin Accessibility Signal Shape us-ing Modified Self-Organizing Maps*. The code in this repository can be used to annotate new chromatin accessibility samples, learn new shapes, or append to a set of existing shapes. We have also included code to replicate our results. Our code is designed for use in a Unix environment and can be run using a command-line interface.
 
## Dependencies
* **Python 3**. Our scripts are all written to run in Python 3. If using an older version of Python, they may be incompatible.
* **bedtools**. This can be installed here: https://github.com/arq5x/bedtools2/releases
* The Unix utilities **shuf**, **cut**, and **awk**. These should be available on most Unix systems.
* The Python modules **numpy**, **scipy**, and **tqdm**.

# Downloading Our Data
To download the data used in our paper, run `common_scripts/download_all_files.sh`. Make sure you run this command from the folder where you wish to download the files. This will download all BAM files, peaks, and ChromHMM mnemonics used in our analyses and concatenate them where needed. This script takes no parameters and will name the BAM files according to cell type.

## Additional Dependencies
* **bamtools**. This can be installed here: https://bioconda.github.io/recipes/bamtools/README.html</li>
* **wget**. This can be installed here: https://ftp.gnu.org/gnu/wget/</li>

# Learning Shapes
The code you will need for this task is in the folder `shape_learning_scripts` and `common_scripts`. Given the input BAM and ChromHMM mnemonic files, follow these steps to learn shapes and associate them with ChromHMM mnemonics.

## Additional Dependencies
* **bamtools**. This can be installed here: https://bioconda.github.io/recipes/bamtools/README.html
* **bigWigToWig**. This can be installed here: https://www.encodeproject.org/software/bigwigtowig/
* **Tensorflow** (for SOM-VN and SOM models). This can be installed here: https://www.tensorflow.org/install/
* **bamCoverage**. This can be installed here: https://deeptools.readthedocs.io/en/develop/content/installation.html
* The python library **pysam**
* **BedOps** (for signal intensity model only). This can be installed here: https://bedops.readthedocs.io/en/latest/index.html

## Steps
1. Download the Kundaje Lab blacklist file from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/ and unzip it using **gunzip**.
2. Run `sort -k1,1 -k2,2n hg38.blacklist.bed > hg38.blacklist_sorted.bed`.
3. Run `bedtools merge -i hg38.blacklist_sorted.bed > hg38.blacklist_merged.bed` to create a merged blacklist.
4. Run `common_scripts/convert_bam_to_wig.sh` to convert the BAM files to chromosome-specific WIG files. This requires the following parameters to be specified:

   * **-b:** The BAM file used as input.
   * **-d:** The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').
   * **-i:** The bin size used to generate the WIG file (default: 10 bp).
   * **-l:** The blacklist regions to exclude. This is the merged blacklist file you created.
   * **-m:** The smoothing length (default: 180 bp).
   * **-n:** The name of the cell line (e.g. Brain). This will be used to locate the input BAM file and to name the output WIG files.
   * **-s:** The file containing a list of chromosome sizes. This is needed for splitting the BAM file by chromosome.   
  
5. Run `common_scripts/create_regions.sh` to create training regions. To create training regions from a permuted WIG file, run `common_scripts/create_regions_permuted.sh`. The requires the following parameters to be specified:

   * **-b:** The bin size used to generate the WIG file (default: 10 bp)
   * **-d:** The directory to contain the split regions
   * **-r:** The region size used for splitting (default: 1 kbp)
   * **-s:** The path to the helper scripts (i.e. the `common_scripts` directory
   * **-w:** The directory containing the WIG file
   
6. Run `shape_learning_scripts/learn_shapes_for_chrom_som_vn.sh` to learn shapes on one chromosome using our method. To learn shapes using other methods, similar scripts are available, with the extensions *_cagt*, *_chromhmmperm*, *_signal*, and *_som*. This requires the following parameters to be specified:

   * **-a:** Directory containing training regions (not needed for *_signal*)
   * **-b:** The bin size used to generate the WIG file (default: 10 bp)
   * **-c:** The chromosome name
   * **-d:** The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/')
   * **-g:** The total number of cells in the SOM grid (only needed for SOM and SOM-VN, default: 100)
   * **-h:** The ChromHMM file used for intersecting
   * **-i:** The number of iterations to use in building the model
   * **-k:** The number of clusters to learn prior to agglomerative clustering (CAGT only, default: 40)
   * **-l:** The learning rate to use in the SOM (default: 0.2)
   * **-m:** The maximum distance for merging to occur in the agglomerative clustering step of CAGT (CAGT only, default: 0.8)
   * **-n:** The single-dimensional neighborhood size to use in the SOM (SOM and SOM-VN only, default: 17)
   * **-p:** The path where the cagt.m file from the CAGT installation is stored (CAGT only)
   * **-r:** The size of the input regions (default: 1000) (not needed for signal model)
   * **-s:** The path to the helper scripts (i.e. the `shape_learning_scripts` directory
   * **-t:** The cutoff to use for cross-correlation significance (not needed for signal model)
   * **-u:** Percentile cutoff file (not needed for signal model)
   * **-w:** The path to the WIG file for this chromosome (signal model only)
   * **-z:** Whether or not you are using PEAS ground truth ("True"/"False")
   
## Description of Helper Scripts
   * `create_index_pysam.py`: Indexes a BAM file using Pysam. This is necessary because, if you index using bamtools, the **bamCoverage** utility cannot process the file.
   * `split_regions.py`: Splits the WIG file into regions to use for training the SOM-VN.
   * `som_vn.py`: The central SOM-VN script. It learnes the shapes given the training regions.
   * `som.py`: The vanilla version of the central SOM-VN script.
   * `permute_chromhmm.py`: Permutes the ChromHMM mnemonics.
   * `permute_wig.py`: Permutes the WIG signal intensities.
   * `extract_signal.py`: Extracts pickled input regions and stores them in CSV files for use by CAGT.
   * run_cagt.m`: Runs CAGT. Note that MATLAB is required to run this script, because CAGT is implemented in MATLAB. CAGT must also be installed; you can download it here: https://github.com/sofiakp/cagt/tree/master/matlab
   * `writeTextResults.m`: This is a hackish solution for running CAGT (There is an I/O issue between our code and the original CAGT). After installing CAGT, you must replace the `writeTextResults.m` file in the `matlab/src` directory with our file. Note that ONLY the I/O format has been changed.
   * `merge_shifted.py`: Consolidates shifted shapes learned by the SOM-VN using cross-correlation.
   * `make_shape_bed.py`: Annotate each region with its closest shape learned by the SOM-VN.
   * `find_chromhmm_distrib.py`: Finds the distribution of ChromHMM mnemonics across each shape.
   * `signal_chromhmm_distrib.py`: Finds the distribution of ChromHMM mnemonics across each signal intensity.

## Output
* A WIG file for each chromosome in the directory `wig` in the base directory.
* Percentile cutoffs for computing crossing count in the `percentile_cutoffs` folder of the directory specified.
* The training regions in the `shifted` directory of the directory specified.
* The shapes learned by the SOM-VN and shifting procedure in the directory `som_vn_output_shifted` in the base directory.
* The training regions annotated with shapes in the directory `som_vn_anno_beds` in the base directory.
* The sorted list of regions annotated with these shapes in the directory `anno_beds_sorted` in the base directory.
* The intersections between our shapes and the ChromHMM regulatory annotations in `som_vn_intersects` in the base directory.
* The shapes with their associated ChromHMM mnemonics in the directory `som_vn_chromhmm_distrib` in the base directory.

# Annotating New Samples
The code you will need for this task is in the folder `annotation_scripts` and `common_scripts`. 

## Steps
1. Run `common_scripts/convert_bam_to_wig.sh` and `common_scripts/create_regions.sh` as described in [Learning Shapes] (#Learning-Shapes).
2. If you want to annotate new samples using region shape, run `python annotation_scripts/make_annotation_bed.py` with the following parameters:

   * The name of the file containing your regions to annotate (generated in step 1).
   * The name of the file containing your learned shapes.
   * The name of the file where you wish to store your annotated BED file.
   * The percentage cutoff for a shape to be associated with a promoter. We use 0.5.
   * The percentage cutoff for a shape to be associated with an enhancer. We use 0.5.
   * The percentage cutoff for a shape to be associated with a repressor. We use. 0.9.
   * The percentage cutoff for a shape to be associated with weak RE. We use 0.9.
   * Whether you want to annotate enhancers only (True / False), as in PEAS.
  
3. If you want to annotate new samples using maximum region signal intensity, run `python annotation_scripts/signal_chromhmm_match.py` with the following parameters:

   * The name of the file containing your regions to annotate (generated in step 1).
   * The name of the file containing your binned signal values.
   * The name of the file where you wish to store your annotated BED file.
   * The percentage cutoff for a shape to be associated with a promoter. We use 0.5.
   * The percentage cutoff for a shape to be associated with an enhancer. We use 0.5.
   * The percentage cutoff for a shape to be associated with a repressor. We use. 0.9.
   * The percentage cutoff for a shape to be associated with weak RE. We use 0.9.
   * Whether you want to annotate enhancers only (True / False), as in PEAS.

# Evaluating Ground Truth on Peaks
The code you will need for this task is in the folder `annotation_scripts`. Here, it is assumed that you have BED files generated using the annotation scripts above, peaks called for your input data, and a BED file with ground truth mnemonics in the format of ChromHMM mnemonics. 

## Steps
1. Run `bedtools intersect -a <annotated-region-bed> -b <peak-bed> > <annotated-peak-bed>`.
2. Run `bedtools intersect -wao -a <annotated-peak-bed> -b <chromhmm-bed> > <peak-chromhmm-overlap-bed>`.
3. Run `python annotation_scripts/save_precision_recall.py <peak-chromhmm-overlap-bed> <precision-recall-file>`.

# Replicating Our Experiments 
Perform the following steps to replicate our experiments:
1. Follow the instructions for [downloading our data](#Downloading-Our-Data).
2. Follow Steps 1-5 from [Learning Shapes](#Learning-Shapes) for each cell type.
3. Run `filter_using_peas.sh` with the following parameters:
   * **-b:** The bin size (default = 10)
   * **-c:** The annotations used as ground truth in the PEAS paper (formatted similarly to ChromHMM).
   * **-r:** The region size (default = 1000)
   * **-t:** The training regions segmented in Step 2.
   * **-w:** The WIG file generated in Step 2.
4. Run `learn_shapes_all_cells_and_chromosomes.sh` with the following parameters:
   * **-c:** The path to the CAGT MATLAB file.
   * **-d:** The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').
   * **-f:** A comma-delimited file with the following format: Cell line name, Training file directory, Permuted training file directory, ChromHMM annotations, Permuted ChromHMM annotations, WIG directory, Training file directory (PEAS), Ground truth  annotations (PEAS), WIG directory (PEAS).
   * **-p:** The project to which you want to charge resources.
   * **-s:** The directory containing the scripts.
5. Run `evaluate_same_chromosome.sh` with the following parameters:
   * **-d:** The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').
   * **-f:** A file listing the inputs for each cell type, which should be in the following format:
        Cell Type, Training Directory, Training Directory (PEAS), ChromHMM Annotations, Ground Truth Annotations from PEAS, Peak File
   * **-p:** The project to which you want to charge resources
   * **-s:** The directory containing the scripts
6. Run `evaluate_cross_chromosome.sh` with the following parameters:
   * **-d:** The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').
   * **-f:** A file listing the inputs for each cell type, which should be in the following format:
        Cell Type, Training Directory, Training Directory (PEAS), ChromHMM Annotations, Ground Truth Annotations from PEAS, Peak File
   * **-p** The project to which you want to charge resources
   * **-s** The directory containing the scripts
7. Run `evaluate_cross_cell_type.sh` with the following parameters:
   * **-d:** The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').
   * **-f:** A file listing the inputs for each cell type, which should be in the following format:
        Cell Type, Training Directory, Training Directory (PEAS), ChromHMM Annotations, Ground Truth Annotations from PEAS, Peak File
   * **-p:** The project to which you want to charge resources
   * **-s:** The directory containing the scripts

# Reproducing Our Figures
To download our data, you will need a system with wget (Most Unix systems should have this). Otherwise, you can download the data manually. You will also need the Python packages glob, pandas, sklearn, matplotlib, and seaborn to run the remaining scripts. Please note that all images will be saved to a file; you do not need a graphical user interface to run this code.
## Fig. 3 
1. Run `python plot_scripts\mnemonic_distrib_per_region.py` on each cell type with the following parameters:
   * The ChromHMM file for that cell type
   * The region size (default: 1000 bp)
   * The output file (CSV format) containing the distribution of each ChromHMM annotation per region
2. Run `python plot_scripts\plot_chromhmm_distribs_density.py` with the following parameters:
   * The CSV file generated in Step 1 for GM12878
   * The CSV file generated in Step 1 for A549
   * The CSV file generated in Step 1 for Brain
   * The CSV file generated in Step 1 for H1
   * The file where you wish to store the figure
   * Whether to remove regions with no annotations (True / False). We removed them.