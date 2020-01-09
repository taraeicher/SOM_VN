#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:10:00 
#!/bin/bash   

#Variables
    USAGE="\n\nThis script is used for associating a set of signals with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift and clustered using k-means. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n\n
    <BIN_SIZE> The bin size used to generate the WIG file (default: 10 bp)\n
    <CHROM> The chromosome name\n
    <BASE_PATH> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <CHROMHMM> The ChromHMM file used for intersecting.\n
    <SCRIPTS> The directory containing the scripts\n
    <WIG> The WIG file for this chromosome\n
    <IS_PEAS> Whether or not this we are using PEAS annotations\n\n"
	
    #Move to the directory containing the scripts.
    cd $SCRIPTS
    
    #Create all needed directories.
    WIGBED="$BASE_PATH/signal_bed"
    if [[ ! -e $WIGBED ]]; then
        mkdir $WIGBED
    fi
    PEAK_INTERSECT_BED="$BASE_PATH/peak_signal_intersects"
    if [[ ! -e $PEAK_INTERSECT_BED ]]; then
        mkdir $PEAK_INTERSECT_BED
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
    bedtools intersect -a $WIGBED/$CHROM.bed -b $PEAKS > $PEAK_INTERSECT_BED/$CHROM.bed
    bedtools intersect -wao -a $PEAK_INTERSECT_BED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python signal_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $CHROMHMM_DISTRIB/$CHROM.pkl $IS_PEAS
        
    #Exit
	wait
	echo Done!
	exit 0