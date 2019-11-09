#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   

#Variables
    USAGE="\n\nThis script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n\n
    <-a> Directory containing training regions\n
    <-b> The bin size used to generate the WIG file (default: 10 bp)\n
    <-c> The chromosome name\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-g> The grid size to use in the SOM (total)\n
    <-h> The ChromHMM file used for intersecting.\n
    <-i> The number of iterations to use in the SOM\n
    <-l> The learning rate to use in the SOM\n
    <-n> The neighborhood size to use in the SOM (single dimensional)\n 
    <-r> The size of the input regions (default: 1000)\n
    <-s> The directory containing the scripts\n
    <-t> The cutoff to use for cross-correlation significance.\n
    <-z> Whether or not this we are using PEAS annotations\n\n"
    
    echo -e $USAGE
    REGION_SIZE=1000
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=10
    TRAINING=""
    CCCUTOFF=0.75
    SCRIPTS=""
    LEARNING_RATE=0.2
    NEIGHBORHOOD=7
    GRID=100
    ITERATIONS=100
    IS_PEAS=False
    while getopts a:b:c:d:g:h:i:l:n:r:s:t:z: option; do
        case "${option}" in
            a) TRAINING=$(realpath $OPTARG);;
            b) BIN_SIZE=$OPTARG;;
            c) CHROM=$OPTARG;;
            d) BASE_PATH=$(realpath $OPTARG);;
            g) GRID=$OPTARG;;
            h) CHROMHMM=$(realpath $OPTARG);;
            i) ITERATIONS=$OPTARG;;
            l) LEARNING_RATE=$OPTARG;;
            n) NEIGHBORHOOD=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            s) SCRIPTS=$(realpath $OPTARG);;
            t) CCCUTOFF=$OPTARG;;
            z) IS_PEAS=$OPTARG;;
        esac
    done
	
    # Move to the directory containing the scripts.
    cd $SCRIPTS
    
    # Create all needed directories.
    SOM_OUT="$BASE_PATH/som_output"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    SOM_SHIFTED="$BASE_PATH/som_som_output_shifted"
    if [[ ! -e $SOM_SHIFTED ]]; then
        mkdir $SOM_SHIFTED
    fi
    ANNOTATED="$BASE_PATH/som_anno_beds"
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
    INTERSECTS="$BASE_PATH/som_intersects"
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    INTERSECTS_SORTED="$BASE_PATH/som_intersects_sorted"
    if [[ ! -e $INTERSECTS_SORTED ]]; then
        mkdir $INTERSECTS_SORTED
    fi
    CHROMHMM_DISTRIB="$BASE_PATH/som_chromhmm_distrib"
    if [[ ! -e $CHROMHMM_DISTRIB ]]; then
        mkdir $CHROMHMM_DISTRIB
    fi

    # Run the SOM.
    python som.py $TRAINING $SOM_OUT/$CHROM.pkl $CUTOFFS $REGION_SIZE $LEARNING_RATE $NEIGHBORHOOD $GRID $ITERATIONS $BIN_SIZE
    echo -e "SOM model is ready for chrom $CHROM.\n"
    
    # Merge shifted regions.
    python merge_shifted.py $SOM_OUT/$CHROM.pkl $SOM_SHIFTED/$CHROM.pkl $CCCUTOFF
    echo -e "Merging complete for chrom $CHROM.\n"
    
    # Annotate regions with shape.
    python make_shape_bed.py $TRAINING $SOM_SHIFTED/$CHROM.pkl $ANNOTATED/$CHROM.bed
    echo -e "Initial annotations complete for chrom $CHROM.\n"
    
    #Intersect regions with ChromHMM.
    bedtools intersect -wao -a $ANNOTATED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $SOM_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl $IS_PEAS
        
    #Exit
	wait
	echo Done!
	exit 0