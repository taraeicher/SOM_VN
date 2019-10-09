#PBS -l nodes=1:ppn=28
#PBS -l walltime=2:00:00 
#!/bin/bash   

#Variables
    USAGE="This script is used for creating training regions from an input WIG file for each chromosome, which can be obtained from a BAM file using convert_bam_to_wig.sh. It splits the WIG file into regions of a specified size with a specified overlap margin and a specified factor for the linear decrease in weights from the center of the region. In this script, a grid of factors and overlap margins are tested, and the optimal pair of parameters is reported.\n
    The following parameters are optional, but recommended:\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-i> The bin size used to generate the WIG file (default: 50 bp)\n
    <-r> The region size used for splitting (default: 4 kbp)\n
    <-w> The directory containing the WIG file\n
    <-o> The directory to contain the split regions\n
    <-m> The margin to use in splitting the regions\n
    <-f> The factor to use in splitting the regions\n
    <-s> The path to the helper scripts"
    
    CELL_LINE=""
    REGION_SIZE=4000
    BASE_FILENAME=""
    CHROMS_NUM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    BIN_SIZE=50
    BLACKLIST=""
    MARGIN=0
    FACTOR=0.5
    WIG=""
    SPLIT_DIR=""
    while getopts n:d:i:r:w:o:m:f:s: option; do
        case "${option}" in
            d) BASE_PATH=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            w) WIG=$OPTARG;;
            o) SPLIT_DIR=$OPTARG;;
            m) MARGIN=$OPTARG;;
            f) FACTOR=$OPTARG;;
            s) SCRIPTS=$(realpath $OPTARG);;
        esac
    done
#Move to location of scripts.
    cd $SCRIPTS
    
#Print message to user.
    echo -e $USAGE
    echo "Creating regions using the following settings:"
    echo "Name is $CELL_LINE"
    echo "Base directory is $BASE_FILENAME"
    echo "Bin size in WIG file is $BIN_SIZE"
    echo "Region size for annotation is $REGION_SIZE"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"
	
#Create all needed directories.
    if [[ ! -e $SPLIT_DIR ]]; then
        mkdir $SPLIT_DIR
    fi

    if [[ ! -e $SPLIT_DIR/shifted/ ]]; then
        mkdir $SPLIT_DIR/shifted/
    fi
    if [[ ! -e $SPLIT_DIR/crossings/ ]]; then
        mkdir $SPLIT_DIR/crossings/
    fi
    if [[ ! -e $SPLIT_DIR/percentile_cutoffs/ ]]; then
    mkdir $SPLIT_DIR/percentile_cutoffs/
    fi
    for c in $CHROMS_NUM;
    do
        python split_regions.py ${WIG}/$c.wig $BIN_SIZE $c $SPLIT_DIR/${c}.pkl $REGION_SIZE $MARGIN $FACTOR 0.95 $SPLIT_DIR/shifted/${c}.pkl $SPLIT_DIR/crossings/${c}.txt $SPLIT_DIR/percentile_cutoffs/${c}.txt
        
        echo "Splitting complete for chromosome $c with factor $factor and margin $margin"

    done    
    
    echo "All splitting complete!"