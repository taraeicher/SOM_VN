#PBS -l nodes=1:ppn=28
#PBS -l walltime=2:00:00 
#!/bin/bash   
# """
# Dependencies:
# 1. Tensorflow, which can be installed here: https://www.tensorflow.org/install/install_linux
# 2. Samtools: https://sourceforge.net/projects/samtools/files/
# 3. Taolib directory: https://github.com/taoliu/taolib
# 4. Gap statistic code: https://github.com/minddrummer/gap
# 5. Added to your path: gap
# 6. Bedtools
# 7. Kent utilities
# 8. Bamtools
# 9. BamCoverage

#Variables
    USAGE="This script is used for creating training regions from an input WIG file for each chromosome, which can be obtained from a BAM file using convert_bam_to_wig.sh. It splits the WIG file into regions of a specified size with a specified overlap margin and a specified factor for the linear decrease in weights from the center of the region. In this script, a grid of factors and overlap margins are tested, and the optimal pair of parameters is reported.\n
    The following parameters are optional, but recommended:\n
    <-b> The bin size used to generate the WIG file (default: 10 bp)\n
    <-c> The count of permutations to run\n
    <-r> The region size used for splitting (default: 1 kbp)\n
    <-w> The directory containing the WIG file\n
    <-d> The directory to contain the split regions"
    
    REGION_SIZE=1000
    CHROMS_NUM="1 2 3" #"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    BIN_SIZE=10
    BLACKLIST=""
    WIG=""
    SPLIT_DIR=""
    COUNT=1
    while getopts b:c:d:r:w: option; do
        case "${option}" in
            b) BIN_SIZE=$OPTARG;;
            c) COUNT=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            w) WIG=$OPTARG;;
            d) SPLIT_DIR=$OPTARG;;
        esac
    done
    
#Print message to user.
    echo -e $USAGE
    echo "Creating regions using the following settings:"
    echo "Base directory is $BASE_FILENAME"
    echo "Bin size in WIG file is $BIN_SIZE"
    echo "Region size for annotation is $REGION_SIZE"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"
	
#Create all needed directories.
    for((i=5;i<=$COUNT;++i)) 
    do
        if [[ ! -e $SPLIT_DIR/$i ]]; then
            mkdir $SPLIT_DIR/$i
        fi

        if [[ ! -e $SPLIT_DIR/$i/shifted/ ]]; then
            mkdir $SPLIT_DIR/$i/shifted/
        fi
        if [[ ! -e $SPLIT_DIR/$i/crossings/ ]]; then
            mkdir $SPLIT_DIR/$i/crossings/
        fi
        if [[ ! -e $SPLIT_DIR/$i/percentile_cutoffs/ ]]; then
        mkdir $SPLIT_DIR/$i/percentile_cutoffs/
        fi
        for c in $CHROMS_NUM;
        do
        
            python permute_wig.py ${WIG}/${c}.wig ${WIG}/${c}_perm_$i.wig
            python split_regions.py ${WIG}/${c}_perm_$i.wig $BIN_SIZE $c $SPLIT_DIR/$i/${c}.pkl $REGION_SIZE 0.95 $SPLIT_DIR/$i/shifted/${c}.pkl $SPLIT_DIR/$i/crossings/${c}.txt $SPLIT_DIR/$i/percentile_cutoffs/${c}.txt
        done
    done    
    
    echo "All splitting complete!"