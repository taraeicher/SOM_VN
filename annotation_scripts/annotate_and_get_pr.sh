#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   

USAGE="\n\nThis script is used for annotating the regions of one cell type and chromosome with learned shapes. Parameters:\n\n
    <-t> The file containing regions to annotate.\n
    <-s> The file containing learned shapes.\n
    <-d> The base directory for output.\n
    <-c> The chromosome number.\n
    <-a> The BED file containing the peaks.\n
    <-h> The ChromHMM BED file.\n
    <-p> The percentage cutoff for associating a shape with a promoter.\n
    <-e> The percentage cutoff for associating a non-promoter shape with an enhancer.\n
    <-r> The percentage cutoff for associating a non-promoter and non-enhancer shape with a repressor.\n
    <-i> The directory containing the scripts\n\n"

echo -e $USAGE
    REGIONS=""
    SHAPES=""
    BASE_PATH=""
    CHROM=""
    CHROMHMM=""
    PROMOTER_CUTOFF=""
    ENHANCER_CUTOFF=""
    REPRESSOR_CUTOFF=""
    SCRIPTS=""
    PEAKS=""
    while getopts t:s:d:c:h:p:e:r:i:a: option; do
        case "${option}" in
            t) REGIONS=$(realpath $OPTARG);;
            s) SHAPES=$(realpath $OPTARG);;
            d) BASE_PATH=$(realpath $OPTARG);;
            c) CHROM=$OPTARG;;
            a) PEAKS=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            p) PROMOTER_CUTOFF=$OPTARG;;
            e) ENHANCER_CUTOFF=$OPTARG;;
            r) REPRESSOR_CUTOFF=$OPTARG;;
            i) SCRIPTS=$(realpath $OPTARG);;
        esac
    done
    
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
    
    # Make the annotated BED file.
    #python make_annotated_bed.py $REGIONS $SHAPES $ANNOTATED_BED/$CHROM.bed $PROMOTER_CUTOFF $ENHANCER_CUTOFF $REPRESSOR_CUTOFF
    
    # Compute the intersect of the annotated regions and the peaks.
    #bedtools intersect -a $ANNOTATED_BED/$CHROM.bed -b $PEAKS > $PEAK_INTERSECT_BED/$CHROM.bed
    bedtools intersect -wao -a $PEAK_INTERSECT_BED/$CHROM.bed -b $CHROMHMM > $CHROMHMM_INTERSECT_BED/$CHROM.bed
    python save_precision_recall.py $CHROMHMM_INTERSECT_BED/$CHROM.bed $PRECISION_RECALL/$CHROM.csv