#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
    #Move to the directory containing the scripts.
    cd $SCRIPTS
    
    #Create all needed directories.
    ANNOTATED_BED="$BASE_PATH/annotated"
    if [[ ! -e $ANNOTATED_BED ]]; then
        mkdir $ANNOTATED_BED
    fi
    PEAK_INTERSECT_BED="$BASE_PATH/peak_intersect"
    if [[ ! -e $PEAK_INTERSECT_BED ]]; then
        mkdir $PEAK_INTERSECT_BED
    fi
    CHROMHMM_INTERSECT_BED="$BASE_PATH/chromhmm_intersect"
    if [[ ! -e $CHROMHMM_INTERSECT_BED ]]; then
        mkdir $CHROMHMM_INTERSECT_BED
    fi
    PRECISION_RECALL="$BASE_PATH/precision_recall"
    if [[ ! -e $PRECISION_RECALL ]]; then
        mkdir $PRECISION_RECALL
    fi
    
    echo $REGIONS
    echo $SHAPES
    echo $ANNOTATED_BED
    echo $CHROM
    
    # Make the annotated BED file.
    python make_annotated_bed.py $REGIONS $SHAPES $ANNOTATED_BED/$CHROM.bed $IS_PEAS
    python ../common_scripts/make_shape_bed.py $REGIONS $SOM_SHIFTED/$CHROM.pkl $ANNOTATED_BED/$CHROM.bed
    
    # Compute the intersect of the annotated regions and the peaks.
    bedtools intersect -a $ANNOTATED_BED/$CHROM.bed -b $PEAKS > $PEAK_INTERSECT_BED/$CHROM.bed
    bedtools intersect -wao -a $PEAK_INTERSECT_BED/$CHROM.bed -b $CHROMHMM > $CHROMHMM_INTERSECT_BED/$CHROM.bed
    python save_precision_recall.py $CHROMHMM_INTERSECT_BED/$CHROM.bed $PRECISION_RECALL/$CHROM.csv $IS_PEAS