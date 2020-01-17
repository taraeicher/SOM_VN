#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00 
#!/bin/bash   

#Variables
    source activate tensorflow_gpu
    module load python/3.6
    USAGE="\n\nThis script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n\n
    <TRAINING> Directory containing training regions\n
    <BIN_SIZE> The bin size used to generate the WIG file (default: 10 bp)\n
    <CHROM> The chromosome name\n
    <BASE_PATH> The base filename where the output files will be stored (e.g. '/root/annoshaperun/').\n
    <GRID> The grid size to use in the SOM (total)\n
    <CHROMHMM> The ChromHMM file used for intersecting.\n
    <ITERATIONS> The number of iterations to use in the SOM\n
    <LEARNING_RATE> The learning rate to use in the SOM\n
    <NEIGHBORHOOD> The neighborhood size to use in the SOM (single dimensional)\n 
    <REGION_SIZE> The size of the input regions (default: 1000)\n
    <SCRIPTS> The directory containing the scripts\n
    <CCCUTOFF> The cutoff to use for cross-correlation significance.\n
    <CUTOFFS> Percentile cutoff file\n
    <IS_PEAS> Whether or not this we are using PEAS annotations\n\n"
    
    # Move to the directory containing the scripts.
    cd $SCRIPTS
    
    # Create all needed directories.
    SOM_OUT="$BASE_PATH/som_vn_output"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    SOM_SHIFTED="$BASE_PATH/som_vn_som_output_shifted"
    if [[ ! -e $SOM_SHIFTED ]]; then
        mkdir $SOM_SHIFTED
    fi
    ANNOTATED="$BASE_PATH/som_vn_anno_beds"
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
    PEAK_INTERSECT_BED="$BASE_PATH/som_vn_peak_intersects"
    if [[ ! -e $PEAK_INTERSECT_BED ]]; then
        mkdir $PEAK_INTERSECT_BED
    fi
    INTERSECTS="$BASE_PATH/som_vn_intersects"
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    INTERSECTS_SORTED="$BASE_PATH/som_vn_intersects_sorted"
    if [[ ! -e $INTERSECTS_SORTED ]]; then
        mkdir $INTERSECTS_SORTED
    fi
    CHROMHMM_DISTRIB="$BASE_PATH/som_vn_chromhmm_distrib"
    if [[ ! -e $CHROMHMM_DISTRIB ]]; then
        mkdir $CHROMHMM_DISTRIB
    fi

    # Run the SOM.
    #python som_vn.py $TRAINING $SOM_OUT/$CHROM.pkl $CUTOFFS $REGION_SIZE $LEARNING_RATE $NEIGHBORHOOD $GRID $ITERATIONS $BIN_SIZE
    #echo -e "SOM model is ready for chrom $CHROM.\n"
    
    # Merge shifted regions.
    #python merge_shifted.py $SOM_OUT/$CHROM.pkl $SOM_SHIFTED/$CHROM.pkl $CCCUTOFF
    #echo -e "Merging complete for chrom $CHROM.\n"
    
    # Annotate regions with shape.
    #python make_shape_bed.py $TRAINING $SOM_SHIFTED/$CHROM.pkl $ANNOTATED/$CHROM.bed
    #echo -e "Initial annotations complete for chrom $CHROM.\n"
    
    # Intersect regions with ChromHMM.
    bedtools intersect -a $ANNOTATED/$CHROM.bed -b $PEAKS > $PEAK_INTERSECT_BED/$CHROM.bed
    bedtools intersect -wao -a $PEAK_INTERSECT_BED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $SOM_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl $IS_PEAS
        
    # Exit
	echo Done!
	exit 0