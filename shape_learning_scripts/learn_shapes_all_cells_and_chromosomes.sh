USAGE="This script is used for annotating each cell type's regions with the regions learned on other cell types. This is done for each chromosome. Parameters:\n
    <-c> The path to the CAGT MATLAB file
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-f> A comma-delimited file containing the following:\n
    \tCell line name,Training file directory,Permuted training file directory,ChromHMM annotations,Permuted ChromHMM annotations,WIG directory,Training file directory (PEAS),Ground truth annotations (PEAS),WIG directory (PEAS)\n
    <-p> The project to which you want to charge resources\n
    <-s> The directory containing the SCRIPTS_DIR\n\n"
    
echo -e $USAGE
BASE_PATH=""
CHROMHMM=""
BIN_SZ=10
REGION_SZ=1000
TRAINING=""
CUTOFF=0.75
SCRIPTS_DIR=""
PARAM_FILE=""
CAGT_P=""
while getopts c:d:f:p:s: option; do
    case "${option}" in
        c) CAGT_P=$(realpath $OPTARG);;
        d) BASE_PATH=$(realpath $OPTARG);;
        f) PARAM_FILE=$(realpath $OPTARG);;
        p) PROJECT=$OPTARG;;
        s) SCRIPTS_DIR=$(realpath $OPTARG);;
    esac
done

# Read in parameters from file.
cell_types=($(awk -F "\"*,\"*" '{print $1;}' $PARAM_FILE | tr '\r\n' ' '))
training_files=($(awk -F "\"*,\"*" '{print $2;}' $PARAM_FILE | tr '\r\n' ' '))
perm_training_files=($(awk -F "\"*,\"*" '{print $3;}' $PARAM_FILE | tr '\r\n' ' '))
chromhmm_anno=($(awk -F "\"*,\"*" '{print $4;}' $PARAM_FILE | tr '\r\n' ' '))
peak=($(awk -F "\"*,\"*" '{print $5;}' $PARAM_FILE | tr '\r\n' ' '))
perm_chromhmm_anno=($(awk -F "\"*,\"*" '{print $6;}' $PARAM_FILE | tr '\r\n' ' '))
wig=($(awk -F "\"*,\"*" '{print $7;}' $PARAM_FILE | tr '\r\n' ' '))
training_files_peas=($(awk -F "\"*,\"*" '{print $8;}' $PARAM_FILE | tr '\r\n' ' '))
ground_truth_peas=($(awk -F "\"*,\"*" '{print $9;}' $PARAM_FILE | tr '\r\n' ' '))
perm_training_files_peas=($(awk -F "\"*,\"*" '{print $10;}' $PARAM_FILE | tr '\r\n' ' '))
perm_ground_truth_peas=($(awk -F "\"*,\"*" '{print $11;}' $PARAM_FILE | tr '\r\n' ' '))
wig_peas=($(awk -F "\"*,\"*" '{print $12;}' $PARAM_FILE | tr '\r\n' ' '))

# Parameters for all algorithms.
grid_size=100
epochs=100
epochs_cagt=1000
cell_type_count=${#cell_types[@]}

# Run for each chromosome
for chrom in 8 #1 8 15 22
do
    # Repeat each run 5 times for robustness.
    for repeat in 1 2 3 4 5
    do       
        # Run for each cell type.
        for (( i=0; i<${cell_type_count}; i++ ));
        do
            # Make training and testing split directories.
            if [[ ! -e ${training_files[$i]}/shifted_split_training_$repeat ]]; then
                mkdir -p ${training_files[$i]}/shifted_split_training_$repeat
            fi
            if [[ ! -e ${training_files[$i]}/percentile_cutoffs_training_$repeat ]]; then
                mkdir -p ${training_files[$i]}/percentile_cutoffs_training_$repeat
            fi
            if [[ ! -e ${training_files[$i]}/shifted_split_testing_$repeat ]]; then
                mkdir -p ${training_files[$i]}/shifted_split_testing_$repeat
            fi
            if [[ ! -e ${training_files[$i]}/percentile_cutoffs_testing_$repeat ]]; then
                mkdir -p ${training_files[$i]}/percentile_cutoffs_testing_$repeat
            fi
            if [[ ! -e ${wig[$i]}_split_training_$repeat ]]; then
                mkdir -p ${wig[$i]}_split_training_$repeat
            fi
            if [[ ! -e ${wig[$i]}_split_testing_$repeat ]]; then
                mkdir -p ${wig[$i]}_split_testing_$repeat
            fi
            if [[ ! -e ${perm_training_files[$i]}/shifted_split_training_$repeat ]]; then
                mkdir -p ${perm_training_files[$i]}/shifted_split_training_$repeat
            fi
            if [[ ! -e ${perm_training_files[$i]}/percentile_cutoffs_training_$repeat ]]; then
                mkdir -p ${perm_training_files[$i]}/percentile_cutoffs_training_$repeat
            fi
            if [[ ! -e ${perm_training_files[$i]}/shifted_split_testing_$repeat ]]; then
                mkdir -p ${perm_training_files[$i]}/shifted_split_testing_$repeat
            fi
            if [[ ! -e ${perm_training_files[$i]}/percentile_cutoffs_testing_$repeat ]]; then
                mkdir -p ${perm_training_files[$i]}/percentile_cutoffs_testing_$repeat
            fi
        
            cutoff=0.95
            # # Randomly choose half of the regions for training.
            # echo ${training_files[$i]} $chrom $repeat
            # python random_split.py ${training_files[$i]}/shifted/${chrom}.pkl ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${training_files[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
            # python ../common_scripts/get_intensity_percentile_regions.py ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl $cutoff ${training_files[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt
            # python ../common_scripts/get_intensity_percentile_regions.py ${training_files[$i]}/shifted_split_testing_$repeat/${chrom}.pkl $cutoff ${training_files[$i]}/percentile_cutoffs_testing_$repeat/${chrom}.txt
            
            # echo ${perm_training_files[$i]} $chrom $repeat
            # python random_split.py ${perm_training_files[$i]}/$repeat/shifted/${chrom}.pkl ${perm_training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${perm_training_files[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
            # python ../common_scripts/get_intensity_percentile_regions.py ${perm_training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl $cutoff ${perm_training_files[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt
            # python ../common_scripts/get_intensity_percentile_regions.py ${perm_training_files[$i]}/shifted_split_testing_$repeat/${chrom}.pkl $cutoff ${perm_training_files[$i]}/percentile_cutoffs_testing_$repeat/${chrom}.txt
            
            
            # python random_split_wig.py ${wig[$i]}/$chrom.wig ${wig[$i]}_split_testing_$repeat/$chrom.wig ${wig[$i]}_split_training_$repeat/$chrom.wig
            
            #if [ ${cell_types[$i]} == "GM12878" ];
            #then
                # if [[ ! -e ${training_files_peas[$i]}/shifted_split_training_$repeat ]]; then
                    # mkdir -p ${training_files_peas[$i]}/shifted_split_training_$repeat
                # fi
                # if [[ ! -e ${training_files_peas[$i]}/percentile_cutoffs_training_$repeat ]]; then
                    # mkdir -p ${training_files_peas[$i]}/percentile_cutoffs_training_$repeat
                # fi
                # if [[ ! -e ${training_files_peas[$i]}/shifted_split_testing_$repeat ]]; then
                    # mkdir -p ${training_files_peas[$i]}/shifted_split_testing_$repeat
                # fi
                # if [[ ! -e ${training_files_peas[$i]}/percentile_cutoffs_testing_$repeat ]]; then
                    # mkdir -p ${training_files_peas[$i]}/percentile_cutoffs_testing_$repeat
                # fi
                # if [[ ! -e ${wig_peas[$i]}_split_training_$repeat ]]; then
                    # mkdir -p ${wig_peas[$i]}_split_training_$repeat
                # fi
                # if [[ ! -e ${wig_peas[$i]}_split_testing_$repeat ]]; then
                    # mkdir -p ${wig_peas[$i]}_split_testing_$repeat
                # fi
                # if [[ ! -e ${perm_training_files_peas[$i]}/shifted_split_training_$repeat ]]; then
                    # mkdir -p ${perm_training_files_peas[$i]}/shifted_split_training_$repeat
                # fi
                # if [[ ! -e ${perm_training_files_peas[$i]}/percentile_cutoffs_training_$repeat ]]; then
                    # mkdir -p ${perm_training_files_peas[$i]}/percentile_cutoffs_training_$repeat
                # fi
                # if [[ ! -e ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat ]]; then
                    # mkdir -p ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat
                # fi
                # if [[ ! -e ${perm_training_files_peas[$i]}/percentile_cutoffs_testing_$repeat ]]; then
                    # mkdir -p ${perm_training_files_peas[$i]}/percentile_cutoffs_testing_$repeat
                # fi
                
                #echo ${training_files_peas[$i]}  $chrom $repeat
                #python random_split.py ${training_files_peas[$i]}/shifted/${chrom}.pkl ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${training_files_peas[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
                #python ../common_scripts/get_intensity_percentile_regions.py ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl $cutoff ${training_files_peas[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt
                #python ../common_scripts/get_intensity_percentile_regions.py ${training_files_peas[$i]}/shifted_split_testing_$repeat/${chrom}.pkl $cutoff ${training_files_peas[$i]}/percentile_cutoffs_testing_$repeat/${chrom}.txt
                
                #echo ${perm_training_files_peas[$i]} $chrom $repeat
                #python random_split.py ${perm_training_files_peas[$i]}/$repeat/shifted/${chrom}.pkl ${perm_training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
                
                #python ../common_scripts/get_intensity_percentile_regions.py ${perm_training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl $cutoff ${perm_training_files_peas[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt
                #python ../common_scripts/get_intensity_percentile_regions.py ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat/${chrom}.pkl $cutoff ${perm_training_files_peas[$i]}/percentile_cutoffs_testing_$repeat/${chrom}.txt
                
                #python random_split_wig.py ${wig_peas[$i]}/$chrom.wig ${wig_peas[$i]}_split_testing_$repeat/$chrom.wig ${wig_peas}_split_training_$repeat/$chrom.wig
            #fi
        
            # Submit all SOM-VN and SOM jobs.
            for alpha_0 in 0.2 0.4 0.6 0.8
            do
                for sigma_0 in 2 4 6 8
                do
                    # # Create a directory.
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_training/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_signalperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_signalperm/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_som/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_som/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    # if [[ ! -e $BASE_PATH/${cell_types[$i]}_som_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        # mkdir -p $BASE_PATH/${cell_types[$i]}_som_training/$chrom/$repeat/$alpha_0/$sigma_0
                    # fi
                    
                    # Normal SOM-VN
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files[$i]}/percentile_cutoffs/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_vn_nogetopts.sh
  
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_vn_nogetopts.sh                    

                    # Permuted WIG SOM-VN
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${perm_training_files[$i]}/$repeat/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_signalperm/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${perm_training_files[$i]}/$repeat/percentile_cutoffs/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_vn_nogetopts.sh
                    
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${perm_training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${perm_training_files[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_vn_nogetopts.sh
                    
                    # Permuted annotation SOM-VN
                    echo $(realpath ${perm_chromhmm_anno[$i]}/${repeat}_chr${chrom}.bed)
                    qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${perm_chromhmm_anno[$i]}/${repeat}_chr${chrom}.bed),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files[$i]}/percentile_cutoffs/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_vn_nogetopts.sh

                    qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${perm_chromhmm_anno[$i]}/${repeat}_chr${chrom}.bed),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_vn_nogetopts.sh
                    
                    # Normal SOM
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_som/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_nogetopts.sh
                    
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_som_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_som_nogetopts.sh
                    
                    # Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} == "GM12878" ];
                    then
                    
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_peas/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_peas/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_peas_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_peas_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_peas_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_peas_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_peas_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_peas_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_somvn_peas_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_somvn_peas_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_som_peas/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_som_peas/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_som_peas_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_som_peas_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
            
                        # Normal SOM-VN with PEAS ground truth
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_peas/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),CUTOFFS=$(realpath ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt),IS_PEAS=True learn_shapes_for_chrom_som_vn_nogetopts.sh
                        
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_peas_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files_peas[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_vn_nogetopts.sh
                        
                        # Permuted WIG SOM-VN with PEAS ground truth
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${perm_training_files_peas[$i]}/$repeat/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${perm_training_files_peas[$i]}/$repeat/percentile_cutoffs/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_vn_nogetopts.sh

                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${perm_training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${perm_training_files_peas[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_vn_nogetopts.sh
                        
                        # Permuted annotation SOM-VN with PEAS ground truth
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_peas_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${perm_ground_truth_peas[$i]}/${repeat}.bed),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_vn_nogetopts.sh
                        
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_somvn_peas_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${perm_ground_truth_peas[$i]}/${repeat}.bed),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,CUTOFFS=$(realpath ${training_files_peas[$i]}/percentile_cutoffs_training_$repeat/${chrom}.txt),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_vn_nogetopts.sh
                        
                        # Normal SOM with PEAS ground truth
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_som_peas/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_nogetopts.sh

                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_som_peas_training/$chrom/$repeat/$alpha_0/$sigma_0/),GRID=$grid_size,CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs,LEARNING_RATE=$alpha_0,NEIGHBORHOOD=$sigma_0,REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_som_nogetopts.sh
                    fi
                done
            done
            
            # Submit all CAGT jobs.
            for k in 20 40 80 160
            do
                for max_dist in 0.3 0.5 0.7 0.9
                do
                    
                    if [[ ! -e $BASE_PATH/${cell_types[$i]}_cagt/$chrom/$repeat/$k/$max_dist ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_cagt/$chrom/$repeat/$k/$max_dist
                    fi
                    if [[ ! -e $BASE_PATH/${cell_types[$i]}_cagt_training/$chrom/$repeat/$k/$max_dist ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_cagt_training/$chrom/$repeat/$k/$max_dist
                    fi
                    
                    # CAGT
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_cagt/$chrom/$repeat/$k/$max_dist/),CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs_cagt,K=$k,MERGE_DIST=$max_dist,CAGT_PATH=$(realpath $CAGT_P),REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_cagt_nogetopts.sh
                    
                    #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_cagt_training/$chrom/$repeat/$k/$max_dist/),CHROMHMM=$(realpath ${chromhmm_anno[$i]}),ITERATIONS=$epochs_cagt,K=$k,MERGE_DIST=$max_dist,CAGT_PATH=$(realpath $CAGT_P),REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_cagt_nogetopts.sh
                    
                    #Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} == "GM12878" ]
                    then
                    
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_cagt_peas/$chrom/$repeat/$k/$max_dist ]]; then
                            mkdir -p $BASE_PATH/${cell_types[$i]}_cagt_peas/$chrom/$repeat/$k/$max_dist
                        fi
                        if [[ ! -e $BASE_PATH/${cell_types[$i]}_cagt_peas_training/$chrom/$repeat/$k/$max_dist ]]; then
                                mkdir -p $BASE_PATH/${cell_types[$i]}_cagt_peas_training/$chrom/$repeat/$k/$max_dist
                        fi
                    
                        #CAGT with PEAS ground truth
                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_cagt_peas/$chrom/$repeat/$k/$max_dist/),CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs_cagt,K=$k,MERGE_DIST=$max_dist,CAGT_PATH=$(realpath $CAGT_P),REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_cagt_nogetopts.sh

                        #qsub -A $PROJECT -m n -v TRAINING=$(realpath ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl),BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_cagt_peas_training/$chrom/$repeat/$k/$max_dist/),CHROMHMM=$(realpath ${ground_truth_peas[$i]}),ITERATIONS=$epochs_cagt,K=$k,MERGE_DIST=$max_dist,CAGT_PATH=$(realpath $CAGT_P),REGION_SIZE=$REGION_SZ,SCRIPTS=$(realpath $SCRIPTS_DIR),CCCUTOFF=$CUTOFF,PEAKS=$(realpath ${peak[$i]}),IS_PEAS=True learn_shapes_for_chrom_cagt_nogetopts.sh
                    fi
                done
            done
            
            # Submit all intensity-only jobs.
            
            if [[ ! -e $BASE_PATH/${cell_types[$i]}_signal/$chrom/$repeat ]]; then
                    mkdir -p $BASE_PATH/${cell_types[$i]}_signal/$chrom/$repeat
            fi
            if [[ ! -e $BASE_PATH/${cell_types[$i]}_signal_training/$chrom/$repeat ]]; then
                    mkdir -p $BASE_PATH/${cell_types[$i]}_signal_training/$chrom/$repeat
            fi
            
            #qsub -A $PROJECT -m n -v BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_signal/$chrom/$repeat/),CHROMHMM=$(realpath ${chromhmm_anno[$i]}),SCRIPTS=$(realpath $SCRIPTS_DIR),WIG=$(realpath $wig/$chrom.wig),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_signal_nogetopts.sh
            
            #qsub -A $PROJECT -m n -v BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_signal_training/$chrom/$repeat/),CHROMHMM=$(realpath ${chromhmm_anno[$i]}),SCRIPTS=$(realpath $SCRIPTS_DIR),WIG=$(realpath ${wig}_split_training_$repeat/$chrom.wig),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_signal_nogetopts.sh
            
            # Submit all intensity-only jobs with PEAS ground truth.
            if [ ${cell_types[$i]} == "GM12878" ]
            then
            
                if [[ ! -e $BASE_PATH/${cell_types[$i]}_peas_signal/$chrom/$repeat ]]; then
                    mkdir -p $BASE_PATH/${cell_types[$i]}_peas_signal/$chrom/$repeat
                fi
                if [[ ! -e $BASE_PATH/${cell_types[$i]}_peas_signal_training/$chrom/$repeat ]]; then
                        mkdir -p $BASE_PATH/${cell_types[$i]}_peas_signal_training/$chrom/$repeat
                fi
            
                #qsub -A $PROJECT -m n -v BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_peas_signal/$chrom/$repeat/),CHROMHMM=$(realpath ${ground_truth_peas[$i]}),SCRIPTS=$(realpath $SCRIPTS_DIR),WIG=$(realpath $wig_peas/$chrom.wig),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_signal_nogetopts.sh
                
                #qsub -A $PROJECT -m n -v BIN_SIZE=$BIN_SZ,CHROM=$chrom,BASE_PATH=$(realpath $BASE_PATH/${cell_types[$i]}_peas_signal_training/$chrom/$repeat/),CHROMHMM=$(realpath ${ground_truth_peas[$i]}),SCRIPTS=$(realpath $SCRIPTS_DIR),WIG=$(realpath ${wig_peas}_split_training_$repeat/$chrom.wig),PEAKS=$(realpath ${peak[$i]}),IS_PEAS=False learn_shapes_for_chrom_signal_nogetopts.sh
            fi
            
        done
    done
done
