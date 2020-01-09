#PBS -l nodes=1:ppn=4
#PBS -l walltime=1:00:00 
#!/bin/bash   

#Variables
    USAGE="\n\nThis script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using CAGT, then merged to correct for signal shift. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n\n
    <TRAINING> Directory containing training regions\n
    <BIN_SIZE> The bin size used to generate the WIG file (default: 10 bp)\n
    <CHROM> The chromosome name\n
    <BASE_PATH> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <CHROMHMM> The ChromHMM file used for intersecting.\n
    <ITERATIONS> The maximum number of iterations for CAGT\n
    <K> The number of clusters to learn prior to agglomerative clustering\n
    <MERGE_DIST> The maximum distance for merging to occur in the agglomerative clustering step of CAGT (default: 0.8)\n
    <CAGT_PATH> The path to the CAGT file\n
    <REGION_SIZE> The size of the input regions (default: 1000)\n
    <SCRIPTS> The directory containing the scripts\n
    <CCCUTOFF> cross-correlation cutoff\n
    <IS_PEAS> Whether or not this we are using PEAS annotations\n\n"

    #Move to the directory containing the scripts.
    cd $SCRIPTS
    
    #Create all needed directories.
    TRAINING_CSV="$BASE_PATH/training_csv"
    if [[ ! -e $TRAINING_CSV ]]; then
        mkdir $TRAINING_CSV
    fi
    MATLAB_MATRIX="$BASE_PATH/cagt_matlab_mat"
    if [[ ! -e $MATLAB_MATRIX ]]; then
        mkdir $MATLAB_MATRIX
    fi
    CAGT_OUT="$BASE_PATH/cagt_output"
    if [[ ! -e $CAGT_OUT ]]; then
        mkdir $CAGT_OUT
    fi
    CAGT_OUT_CSV="$BASE_PATH/cagt_output_csv"
    if [[ ! -e $CAGT_OUT_CSV ]]; then
        mkdir $CAGT_OUT_CSV
    fi
    CAGT_SHIFTED="$BASE_PATH/cagt_output_shifted"
    if [[ ! -e $CAGT_SHIFTED ]]; then
        mkdir $CAGT_SHIFTED
    fi
    ANNOTATED="$BASE_PATH/cagt_anno_beds"
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
    PEAK_INTERSECT_BED="$BASE_PATH/cagt_peak_intersects"
    if [[ ! -e $PEAK_INTERSECT_BED ]]; then
        mkdir $PEAK_INTERSECT_BED
    fi
    INTERSECTS="$BASE_PATH/cagt_intersects"
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    INTERSECTS_SORTED="$BASE_PATH/cagt_intersects_sorted"
    if [[ ! -e $INTERSECTS_SORTED ]]; then
        mkdir $INTERSECTS_SORTED
    fi
    CHROMHMM_DISTRIB="$BASE_PATH/cagt_chromhmm_distrib"
    if [[ ! -e $CHROMHMM_DISTRIB ]]; then
        mkdir $CHROMHMM_DISTRIB
    fi

    #Extract the signal and run CAGT.
    python extract_signal.py $TRAINING $TRAINING_CSV/$CHROM.csv
    module load matlab
    matlab -nodisplay -nodesktop -r "run_cagt('$TRAINING_CSV/$CHROM.csv','$MATLAB_MATRIX','$CHROM','$CAGT_OUT_CSV/$CHROM.csv','$CAGT_PATH', '$MERGE_DIST', '$ITERATIONS', '$K')"
    python convert_to_pickle.py $CAGT_OUT_CSV/$CHROM.csv $CAGT_OUT/$CHROM.pkl
    echo -e "CAGT model is ready for chrom $CHROM.\n"
    
    #Merge shifted regions.
    python merge_shifted.py $CAGT_OUT/$CHROM.pkl $CAGT_SHIFTED/$CHROM.pkl $CCCUTOFF
    echo -e "Merging complete for chrom $CHROM.\n"
    
    #Annotate regions with shape.
    python make_shape_bed.py $TRAINING $CAGT_SHIFTED/$CHROM.pkl $ANNOTATED/$CHROM.bed
    echo -e "Initial annotations complete for chrom $CHROM.\n"
    
    #Intersect regions with ChromHMM.
    bedtools intersect -a $ANNOTATED/$CHROM.bed -b $PEAKS > $PEAK_INTERSECT_BED/$CHROM.bed
    bedtools intersect -wao -a $PEAK_INTERSECT_BED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $CAGT_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl $IS_PEAS
        
    #Exit
	wait
	echo Done!
	exit 0