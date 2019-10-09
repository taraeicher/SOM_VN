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
    <p>The code you will need for this task is in the folder <b>shape_learning_scripts</b> and <b>common_scripts</b>. Given the input BAM and ChromHMM mnemonic files, follow these steps to learn shapes and associate them with ChromHMM mnemonics.</p>
    <h3>Additional Dependencies</h3>
        <ul>
            <li><b>bamtools</b>. This can be installed here: https://bioconda.github.io/recipes/bamtools/README.html</li>
            <li><b>bigWigToWig</b>. This can be installed here: https://www.encodeproject.org/software/bigwigtowig/</li>
            <li><b>Tensorflow</b> (for VNSSOM and Vanilla SOM models). This can be installed here: https://www.tensorflow.org/install/.</li>
            <li><b>bamCoverage</b>. This can be installed here: https://deeptools.readthedocs.io/en/develop/content/installation.html</li>
            <li>The python library pysam</li>
            <li><b>BedOps</b> (for signal intensity model only). This can be installed here: https://bedops.readthedocs.io/en/latest/index.html</li>
        </ul>
    <h3>Steps</h3>
        <ol>
            <li>Download the Kundaje Lab blacklist file from http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/</li>
            <li>Run <b>bedtools merge -i hg38.blacklist.bed -o hg38.blacklist_merged.bed</b> to create a merged blacklist.</li>
            <li>Run <b>shape_learning_scripts/convert_bam_to_wig.sh</b> to convert the BAM files to chromosome-specific WIG files. This requires the following parameters to be specified:
                <ul>
                    <li> <b>-n:</b> The name of the cell line (e.g. Brain). This will be used to locate the input BAM file and to name the output WIG files.</li>
                    <li> <b>-d:</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
                    <li> <b>-b:</b> The BAM file used as input.</li>
                    <li> <b>-i:</b> The bin size used to generate the WIG file (default: 50 bp).</li>
                    <li> <b>-s:</b> The file containing a list of chromosome sizes. This is needed for splitting the BAM file by chromosome.</li>
                    <li> <b>-l:</b> The blacklist regions to exclude. This is the merged blacklist file you created.</li>
                </ul>
            </li>
            <li> Run <b>common_scripts/create_regions.sh</b> to create training regions. To create training regions from a permuted WIG file, run <b>common_scripts/create_regions_permuted.sh</b>. The requires the following parameters to be specified:
                <ul>
                    <li><b>-n</b> The name of the cell line (e.g. Brain)</li>
                    <li><b>-d</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
                    <li><b>-i</b> The bin size used to generate the WIG file (default: 50 bp)</li>
                    <li><b>-r</b> The region size used for splitting (default: 4 kbp)</li>
                    <li><b>-w</b> The directory containing the WIG file</li>
                    <li><b>-o</b> The directory to contain the split regions</li>
                    <li><b>-m</b> The margin to use in splitting the regions</li>
                    <li><b>-f</b> The factor to use in splitting the regions</li>
                    <li><b>-s</b> The path to the helper scripts (i.e. the <b>common_scripts</b> directory.</li>
                </ul>
            </li>
            <li> Run <b>shape_learning_scripts/learn_shapes_for_chrom_vnssom.sh</b> to learn shapes on one chromosome using our method. To learn shapes using other methods, similar scripts are available, with the extensions <b>_cagt</b>, <b>_chromhmmperm</b>, <b>_signal</b>, and <b>_som</b>. This requires the following parameters to be specified:
                <ul>
                    <li><b>-d:</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
                    <li><b>-h:</b> The ChromHMM file used for intersecting.</li>
                    <li><b>-i:</b> The bin size used to generate the WIG file (default: 50 bp)</li>
                    <li><b>-r:</b> The size of the input regions (default: 4000) (not needed for <b>_signal</b>)</li>
                    <li><b>-t:</b> The cutoff to use for cross-correlation significance (not needed for <b>_signal</b>)</li>
                    <li><b>-a:</b> Directory containing training regions (not needed for <b>_signal</b>)</li>
                    <li><b>-u:</b> Percentile cutoff file (not needed for <b>_signal</b>)</li>
                    <li><b>-c:</b> The chromosome name</li>
                    <li><b>-p:</b> The path where the cagt.m file from the CAGT installation is stored (<b>_cagt</b> only)</li>
                    <li><b>-w:</b> The path to the WIG file for this chromosome (<b>_signal</b> only).</li>
                    <li><b>-s:</b> The path to the helper scripts (i.e. the <b>shape_learning_scripts</b> directory.</li>
                </ul>
            Note that you may also learn shapes on all 24 chromosomes in parallel using the scripts <b>shape_learning_scripts/run_all_chroms_vnssom.sh</b>, <b>_som.sh</b>, <b>_signal.sh</b>, or <b>_cagt.sh</b> if you have access to a supercomputer that uses the PBS queueing system. This will run using the default parameters in our paper. In this case, the only parameters you need are:
                <ul>
                    <li><b>-d:</b> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').</li>
                    <li><b>-h:</b> The ChromHMM file used for intersecting.</li>
                    <li><b>-a:</b> Directory containing training regions (not needed for <b>_signal</b>)</li>
                    <li><b>-c:</b> The chromosome name</li>
                    <li><b>-p:</b> The project number to which you wish to charge your run</li>
                    <li><b>-g:</b> The path where the cagt.m file from the CAGT installation is stored (<b>_cagt</b> only)</li>
                    <li><b>-w:</b> The path to the WIG file for this chromosome (<b>_signal</b> only).</li>
                    <li><b>-s:</b> The path to the helper scripts (i.e. the <b>shape_learning_scripts</b> directory.</li>
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
            <li><b>writeTextResults.m:</b> This is a hackish solution for running CAGT (There is a formatting issue with CAGT that cannot be resolved at this time). After installing CAGT, you must replace the <b>writeTextResults.m</b> file in the <b>matlab/src</b> directory with our file.</li>
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
<h1>Annotating New Samples</h1>
    <p>The code you will need for this task is in the folder <b>annotation_scripts</b> and <b>common_scripts</b>. 
    <h3>Steps</h3>
        <ol>
            <li> Run <b>common_scripts/convert_bam_to_wig.sh</b> and <b>common_scripts/create_regions.sh</b> as described in the <b>Learning Shapes</b> section.</li>
            <li>If you want to annotate new samples using region shape, run <b>python annotation_scripts/make_annotation_bed.py</b> with the following parameters:
                <ul>
                    <li>The name of the file containing your regions to annotate (generated in step 1).</li>
                    <li>The name of the file containing your learned shapes (as provided in <b>learned_shapes</b>).</li>
                    <li>The name of the file where you wish to store your annotated BED file.</li>
                    <li>The percentage cutoff for a shape to be associated with a promoter</li>
                    <li>The percentage cutoff for a shape to be associated with an enhancer, given that it is not a promoter</li>
                    <li>The percentage cutoff for a shape to be associated with a repressor, given that it is neither a promoter nor an enhancer</li>
                </ul>
            </li>
            <li>If you want to annotate new samples using maximum region signal intensity, run <b>python annotation_scripts/signal_chromhmm_match.py</b> with the following parameters:
                <ul>
                    <li>The name of the file containing your regions to annotate (generated in step 1).</li>
                    <li>The name of the file containing your binned signal values (as provided in <b>learned_shapes</b>).</li>
                    <li>The name of the file where you wish to store your annotated BED file.</li>
                    <li>The percentage cutoff for a region's maximum intensity to be associated with a promoter</li>
                    <li>The percentage cutoff for a region's maximum intensity to be associated with an enhancer, given that it is not a promoter</li>
                    <li>The percentage cutoff for a region's maximum intensity to be associated with a repressor, given that it is neither a promoter nor an enhancer</li>
                </ul>
            </li>
        </ol>
<h1>Evaluating Ground Truth on Peaks</h1>
    <p>The code you will need for this task is in the folder <b>annotation_scripts</b>. Here, it is assumed that you have BED files generated using the annotation scripts above, peaks called for your input data, and a BED file with ground truth mnemonics in the format of ChromHMM mnemonics. 
    <h3>Steps</h3>
        <ol>
            <li> Run <b>bedtools intersect -a <i>annotated-region-bed</i> -b <i>peak-bed</i> > <i>annotated-peak-bed</i></b></li>
            <li> Run <b>bedtools intersect -wao -a <i>annotated-peak-bed</i> -b <i>chromhmm-bed</i> > <i>peak-chromhmm-overlap-bed</i></b></li>
            <li>Run <b>python annotation_scripts/save_precision_recall.py <i>peak-chromhmm-overlap-bed</i> <i>precision-recall-file</i></b></li>
        </ol>
<h1>Reproducing Our Figures</h1>
    <p>To download our data, you will need a system with wget (Most Unix systems should have this). Otherwise, you can download the data manually. You will also need the Python packages glob, pandas, sklearn, matplotlib, and seaborn to run the remaining scripts. Please note that all images will be saved to a file; you do not need a graphical user interface to run this code.</p>
    <ul>
        <li><b>Fig. 3: </b> Run <b>plot_precision_recall.py</b> followed by <b>plot_precision_recall_all.py</b> (for part A) and <b>plot_true_distribs_all.py</b> (for part B).</li>
        <li><b>Fig. 4: </b> Run <b>run_all_chromosome_iterations.sh</b> on each cell type to generate the density plot followed by <b>save_precision_recall.py</b> and <b>plot_precision_recall_densities.py</b>.</li>
        <li><b>Fig. 5: </b> Run <b>plot_precision_recall_nobaselines.py</b>.</li>
        <li><b>Fig. 6: </b> Run <b>plot_precision_recall_nopromoter_abovethreshonly.py</b>.</li>
        <li><b>Supplementary Fig. 1: </b> Run <b>plot_wig_distribs_violin.py</b>.</li>
        <li><b>Supplementary Fig. 4: </b> Run <b>plot_chromhmm_distribs_violin.py</b>.</li>
        <li><b>Supplementary Fig. 5: </b> Run <b>annotation_similarity_heatmap.py</b>.</li>
        <li><b>Supplementary Fig. 6: </b> Run <b>print_annotated_shapes.py</b>.</li>
        <li><b>Supplementary Fig. 7, 8, and 9: </b> Run <b>plot_precision_recall.py</b>.</li>
        <li><b>Supplementary Fig. 10: </b> Run <b>plot_crosscorr_distrib.py</b>.</li>
        <li><b>Supplementary Fig. 11: </b> Run <b>plot_precision_recall_nobaselines.py</b>.</li>
    </ul>