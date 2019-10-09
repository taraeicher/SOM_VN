USAGE="This script is used for associating signal intensity with an RE using ChromHMM. This script will run the training on all 24 chromosomes. Parameters:\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-w> WIG file.\n
    <-c> The chromosome name\n
    <-p> The project to which you want to charge resources\n
    <-s> The directory containing the scripts"
    
    echo -e $USAGE
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=50
    TRAINING=""
    SCRIPTS=""
    while getopts h:d:c:a:p:w:s: option; do
        case "${option}" in
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            c) CHROM=$(realpath $OPTARG);;
            p) PROJECT=$(realpath $OPTARG);;
            s) SCRIPTS=$(realpath $OPTARG);;
        esac
    done

for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -d $BASE_PATH -h $CHROMHMM -c $chrom -u $BASE_PATH/percentilecutoffs/$chrom.txt -s $SCRIPTS
done
