#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   

#Variables
    USAGE="\n\nThis script is used for associating a set of signals with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift and clustered using k-means. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n\n
    <-b> The bin size used to generate the WIG file (default: 10 bp)\n
    <-c> The chromosome name\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-s> The directory containing the scripts\n
    <-w> The WIG file for this chromosome\n
    <-z> Whether or not this we are using PEAS annotations\n\n"
    
    echo -e $USAGE
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=10
    CCCUTOFF=0.75
    SCRIPTS=""
    IS_PEAS=False
    while getopts b:c:d:h:s:w:z: option; do
        case "${option}" in
            b) BIN_SIZE=$OPTARG;;
            c) CHROM=$OPTARG;;
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            s) SCRIPTS=$(realpath $OPTARG);;
            w) WIG=$(realpath $OPTARG);;
            z) IS_PEAS=$OPTARG;;
        esac
    done
	
    #Move to the directory containing the scripts.
    cd $SCRIPTS
    
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
    wig2bed --zero-indexed < $WIG > $WIGBED/$CHROM.bed
    
    # Intersect signal BED file with ChromHMM.
    bedtools intersect -wao -a $WIGBED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python signal_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $CHROMHMM_DISTRIB/$CHROM.pkl $IS_PEAS
        
    #Exit
	wait
	echo Done!
	exit 0