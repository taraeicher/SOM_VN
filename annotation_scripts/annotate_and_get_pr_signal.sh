#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   

USAGE="\n\nThis script is used for annotating the regions of one cell type and chromosome with learned shapes. Parameters:\n\n
    <-a> The BED file containing the peaks.\n
    <-c> The chromosome number.\n
    <-d> The base directory for output.\n
    <-e> The percentage cutoff for associating a non-promoter shape with an enhancer.\n
    <-h> The ChromHMM BED file.\n
    <-i> The directory containing the scripts\n
    <-p> The percentage cutoff for associating a shape with a promoter.\n
    <-r> The percentage cutoff for associating a non-promoter and non-enhancer shape with a repressor.\n
    <-s> The file containing signal intensity associations.\n
    <-t> The file containing regions to annotate.\n
    <-z> Whether or not this we are using PEAS annotations\n"

echo -e $USAGE
    REGIONS=""
    SIGNALS=""
    BASE_PATH=""
    CHROM=""
    CHROMHMM=""
    PROMOTER_CUTOFF=""
    ENHANCER_CUTOFF=""
    REPRESSOR_CUTOFF=""
    SCRIPTS=""
    PEAKS=""
    IS_PEAS=False
    while getopts a:c:d:e:h:i:p:r:s:t:w:z: option; do
        case "${option}" in
            a) PEAKS=$(realpath $OPTARG);;
            c) CHROM=$OPTARG;;
            d) BASE_PATH=$(realpath $OPTARG);;
            e) ENHANCER_CUTOFF=$OPTARG;;
            h) CHROMHMM=$(realpath $OPTARG);;
            i) SCRIPTS=$(realpath $OPTARG);;
            p) PROMOTER_CUTOFF=$OPTARG;;
            r) REPRESSOR_CUTOFF=$OPTARG;;
            s) SIGNALS=$(realpath $OPTARG);;
            t) REGIONS=$(realpath $OPTARG);;
            w) WEAK_CUTOFF=$OPTARG;;
            z) IS_PEAS=$OPTARG;;
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
    python signal_chromhmm_match.py $REGIONS $SIGNALS $ANNOTATED_BED/$CHROM.bed $PROMOTER_CUTOFF $ENHANCER_CUTOFF $REPRESSOR_CUTOFF $WEAK_CUTOFF $IS_PEAS
    
    # Compute the intersect of the annotated regions and the peaks.
    bedtools intersect -a $ANNOTATED_BED/$CHROM.bed -b $PEAKS > $PEAK_INTERSECT_BED/$CHROM.bed
    bedtools intersect -wao -a $PEAK_INTERSECT_BED/$CHROM.bed -b $CHROMHMM > $CHROMHMM_INTERSECT_BED/$CHROM.bed
    python save_precision_recall.py $CHROMHMM_INTERSECT_BED/$CHROM.bed $PRECISION_RECALL/$CHROM.csv $IS_PEAS