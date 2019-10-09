#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   

#Variables
    USAGE="This script is used for associating a set of signals with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift and clustered using k-means. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-i> The bin size used to generate the WIG file (default: 50 bp)\n
    <-c> The chromosome name
    <-w> The WIG directory for this chromosome"
    
    echo -e $USAGE
    REGION_SIZE=4000
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=50
    CCCUTOFF=0.75
    while getopts h:d:c:i:w: option; do
        case "${option}" in
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            w) WIG=$(realpath $OPTARG);;
            c) CHROM=$(realpath $OPTARG);;
        esac
    done
	
    for 
    
    #Create all needed directories.
    WIGBED="$BASE_PATH/signal_bed"
    if [[ ! -e $WIGBED ]]; then
        mkdir $WIGBED
    fi
    INTERSECTS="$BASE_PATH/signal_intersects"
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    INTERSECTS_SORTED="$BASE_PATH/signal_intersects_sorted"
    if [[ ! -e $INTERSECTS_SORTED ]]; then
        mkdir $INTERSECTS_SORTED
    fi
    CHROMHMM_DISTRIB="$BASE_PATH/signal_chromhmm_distrib"
    if [[ ! -e $CHROMHMM_DISTRIB ]]; then
        mkdir $CHROMHMM_DISTRIB
    fi
    
    # Convert WIG files to BED files.
    wig2bed < $WIG/$CHROM.wig > $WIGBED/$CHROM.bed
    
    # Intersect signal BED file with ChromHMM.
    bedtools intersect -wao -a $WIGBED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $SOM_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl
        
    #Exit
	wait
	echo Done!
	exit 0