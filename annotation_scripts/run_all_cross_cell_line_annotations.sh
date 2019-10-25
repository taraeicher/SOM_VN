USAGE="This script is used for annotating each cell type's regions with the regions learned on other cell types. This is done for each chromosome. Parameters:\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-t> The cutoff to use for cross-correlation significance.\n
    <-a> Directory containing training regions\n
    <-c> The chromosome name\n
    <-p> The project to which you want to charge resources\n
    <-s> The directory containing the scripts"
    
    echo -e $USAGE
    BASE_PATH=""
    BASE_OUTPUT_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=10
    TRAINING=""
    CCCUTOFF=0.75
    SCRIPTS=""
    while getopts h:d:c:a:p:s:o: option; do
        case "${option}" in
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            a) TRAINING=$(realpath $OPTARG);;
            c) CHROM=$(realpath $OPTARG);;
            p) PROJECT=$(realpath $OPTARG);;
            s) SCRIPTS=$(realpath $OPTARG);;
            o) BASE_OUTPUT_PATH=$(realpath $OPTARG);;
        esac
    done

for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    for cell_src in "GM12878" "Brain" "A549" "H1"
    do
        for cell_target in "GM12878" "Brain" "A549" "H1"
        do
            for p_promoter in 1 2 3 4 5 6 7 8 9
            do
                for p_enhancer in 1 2 3 4 5 6 7 8 9
                do
                    for p_repressor in 1 2 3 4 5 6 7 8 9
                    do
                        qsub -A $PROJECT annotate_and_get_pr.sh -t $BASE_PATH/$cell_target/$TRAINING/$chrom.pkl -s $BASE_PATH/$cell_src/$SHAPES/$chrom.pkl -d $BASE_OUTPUT_PATH/${cell_src}_${cell_target}_${p_promoter}_${p_enhancer}_${p_repressor} -c $chrom -h $CHROMHMM/$cell_target.bed -p 0.$p_promoter -e 0.$p_enhancer -r 0.$p_repressor -i $SCRIPTS                        
                    done
                done
            done
        done
    done
done
