USAGE="This script is used for annotating each cell type's regions with the regions learned on other cell types. This is done for each chromosome. Parameters:\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-f> A file listing the inputs for each cell type, which should be in the following format:\n
       Cell Type, Training Directory, Training Directory (PEAS), ChromHMM Annotations, Ground Truth Annotations from PEAS, Peak File\n
    <-p> The project to which you want to charge resources\n
    <-s> The directory containing the scripts"
    
echo -e $USAGE
BASE_PATH=""
SCRIPTS=""
PARAM_FILE=""
PROJECT=""
while getopts d:f:p:s: option; do
    case "${option}" in
        d) BASE_PATH=$(realpath $OPTARG);;
        f) PARAM_FILE=$(realpath $OPTARG);;
        p) PROJECT=$OPTARG;;
        s) SCRIPTS=$(realpath $OPTARG);;
    esac
done
    
# Read in parameters from file.
cell_types=$(awk -F "\"*,\"*" '{print $1;}' $PARAM_FILE | tr '\r\n' ' ')
training_files=$(awk -F "\"*,\"*" '{print $2;}' $PARAM_FILE | tr '\r\n' ' ')
training_files_peas=$(awk -F "\"*,\"*" '{print $3;}' $PARAM_FILE | tr '\r\n' ' ')
chromhmm_anno=$(awk -F "\"*,\"*" '{print $4;}' $PARAM_FILE | tr '\r\n' ' ')
chromhmm_anno_peas=$(awk -F "\"*,\"*" '{print $5;}' $PARAM_FILE | tr '\r\n' ' ')
peaks=$(awk -F "\"*,\"*" '{print $6;}' $PARAM_FILE | tr '\r\n' ' ')

cell_type_count=${#cell_types[@]}
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    for (( i=1; i<${cell_type_count}+1; i++ ));
    do
        for repeat in 1 2 3 4 5
        do
            for alpha_0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
            do
                for sigma_0 in 1 2 3 4 5 6 7 8 9
                do
                    qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_vn/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_vn_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z False

                   qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_vn_signalperm/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_vn_signalperm_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z False 
                   
                   qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_vn_chromhmmperm/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_vn_chromhmmperm_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS  -z False
                   
                   qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z False
                   
                   if [ ${cell_types[$i]} = "GM12878"]
                    then
                        qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_vn_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_vn_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno_peas[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z True

                        qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_vn_signalperm_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_vn_signalperm_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno_peas[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z True
                   
                        qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_vn_chromhmmperm_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_vn_chromhmmperm_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno_peas[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS  -z True
                   
                        qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/som_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/ -d $BASE_PATH/som_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -c $chrom -h ${chromhmm_anno_peas[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z True
                   fi
                done
            done
            # Submit all CAGT jobs.
            for k in 20 40 80 160
            do
                for max_dist in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
                do
                    # CAGT
                    qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/cagt/${cell_types[$i]}_split_training/$chrom/$repeat/$k/$max_dist/ -d $BASE_PATH/cagt_intrachrom/${cell_types[$i]}/$chrom/$repeat/$k/$max_dist/ -c $chrom -h ${chromhmm_anno[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z False

                    # Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} = "GM12878"]
                    then
                        # CAGT with PEAS ground truth
                        qsub -A $PROJECT annotate_and_get_pr.sh -a ${peaks[$i]} -t ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/cagt_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$k/$max_dist/ -d $BASE_PATH/cagt_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$k/$max_dist/ -c $chrom -h ${chromhmm_anno_peas[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z True
                    fi
                done
            done
            
            qsub -A $PROJECT annotate_and_get_pr_signal.sh -a ${peaks[$i]} -t ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/signal/${cell_types[$i]}_split_training/$chrom/ -d $BASE_PATH/signal_intrachrom/${cell_types[$i]}/$chrom/ -c $chrom -h ${chromhmm_anno[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z False
            
            if [ ${cell_types[$i]} = "GM12878"]
            then
                qsub -A $PROJECT annotate_and_get_pr_signal.sh -a ${peaks[$i]} -t ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl -s $BASE_PATH/signal_peas/${cell_types[$i]}_split_training/$chrom/ -d $BASE_PATH/signal_peas_intrachrom/${cell_types[$i]}/$chrom/ -c $chrom -h ${chromhmm_anno_peas[$i]} -p 0.5 -e 0.5 -r 0.9 -w 0.9 -i $SCRIPTS -z False
            fi
        done
    done
done
