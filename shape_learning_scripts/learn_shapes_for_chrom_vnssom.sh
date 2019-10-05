#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   

#Variables
    USAGE="This script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift and clustered using k-means. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-i> The bin size used to generate the WIG file (default: 50 bp)\n
    <-r> The size of the input regions (default: 4000)\n
    <-t> The cutoff to use for cross-correlation significance.\n
    <-a> Directory containing training regions\n
    <-u> Percentile cutoff file
    <-c> The chromosome name"
    
    echo -e $USAGE
    REGION_SIZE=4000
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=50
    TRAINING=""
    CCCUTOFF=0.75
    while getopts h:d:c:i:r:t:a:u: option; do
        case "${option}" in
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            a) TRAINING=$(realpath $OPTARG);;
            u) CUTOFFS=$(realpath $OPTARG);;
            c) CHROM=$(realpath $OPTARG);;
            t) CCCUTOFF=$OPTARG;;
        esac
    done
	
    for 
    
    #Create all needed directories.
    SOM_OUT="$BASE_PATH/vnssom_output"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    SOM_SHIFTED="$BASE_PATH/vnssom_som_output_shifted"
    if [[ ! -e $SOM_SHIFTED ]]; then
        mkdir $SOM_SHIFTED
    fi
    ANNOTATED="$BASE_PATH/vnssom_anno_beds"
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
    INTERSECTS="$BASE_PATH/vnssom_intersects"
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    INTERSECTS_SORTED="$BASE_PATH/vnssom_intersects_sorted"
    if [[ ! -e $INTERSECTS_SORTED ]]; then
        mkdir $INTERSECTS_SORTED
    fi
    CHROMHMM_DISTRIB="$BASE_PATH/vnssom_chromhmm_distrib"
    if [[ ! -e $CHROMHMM_DISTRIB ]]; then
        mkdir $CHROMHMM_DISTRIB
    fi

    #Run the SOM.
    python vnssom.py $TRAINING $SOM_OUT/$CHROM.pkl $CUTOFFS $REGION_SIZE $BIN_SIZE
    echo -e "SOM model is ready for chrom $CHROM.\n"
    
    #Merge shifted regions.
    python merge_shifted.py $SOM_OUT/$CHROM.pkl $SOM_SHIFTED/$CHROM.pkl $CCCUTOFF
    echo -e "Merging complete for chrom $CHROM.\n"
    
    #Annotate regions with shape.
    python make_shape_bed.py $TRAINING/$CHROM.pkl $SOM_SHIFTED/$CHROM.pkl $ANNOTATED/$CHROM.pkl
    echo -e "Initial annotations complete for chrom $CHROM.\n"
    
    #Intersect regions with ChromHMM.
    bedtools intersect -wao -a $ANNOTATED/$CHROM.pkl -b $CHROMHMM > $INTERSECTS/$CHROM.pkl
    bedtools sort -i $INTERSECTS/$CHROM.pkl > $INTERSECTS_SORTED/$CHROM.pkl
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.pkl $SOM_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl
        
    #Exit
	wait
	echo Done!
	exit 0