USAGE="This script is used for annotating each cell type's regions with the regions learned on other cell types. This is done for each chromosome. Parameters:\n
    <-c> The path to the CAGT MATLAB file
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-f> A comma-delimited file containing the following:\n
    \tCell line name,Training file directory,Permuted training file directory,ChromHMM annotations,Permuted ChromHMM annotations,WIG directory,Training file directory (PEAS),Ground truth annotations (PEAS),WIG directory (PEAS)\n
    <-p> The project to which you want to charge resources\n
    <-s> The directory containing the scripts\n\n"
    
echo -e $USAGE
BASE_PATH=""
CHROMHMM=""
BIN_SIZE=10
REGION_SIZE=1000
TRAINING=""
CCCUTOFF=0.75
SCRIPTS=""
PARAM_FILE=""
CAGT_PATH=""
while getopts c:d:f:p:s: option; do
    case "${option}" in
        c) CAGT_PATH=$(realpath $OPTARG);;
        d) BASE_PATH=$(realpath $OPTARG);;
        f) PARAM_FILE=$(realpath $OPTARG);;
        p) PROJECT=$(realpath $OPTARG);;
        s) SCRIPTS=$(realpath $OPTARG);;
    esac
done

# Read in parameters from file.
cell_types=$(awk -F "\"*,\"*" '{print $1;}' $PARAM_FILE | tr '\r\n' ' ')
training_files=$(awk -F "\"*,\"*" '{print $2;}' $PARAM_FILE | tr '\r\n' ' ')
perm_training_files=$(awk -F "\"*,\"*" '{print $3;}' $PARAM_FILE | tr '\r\n' ' ')
chromhmm_anno=$(awk -F "\"*,\"*" '{print $4;}' $PARAM_FILE | tr '\r\n' ' ')
perm_chromhmm_anno=$(awk -F "\"*,\"*" '{print $5;}' $PARAM_FILE | tr '\r\n' ' ')
wig=$(awk -F "\"*,\"*" '{print $6;}' $PARAM_FILE | tr '\r\n' ' ')
training_files_peas=$(awk -F "\"*,\"*" '{print $7;}' $PARAM_FILE | tr '\r\n' ' ')
ground_truth_peas=$(awk -F "\"*,\"*" '{print $8;}' $PARAM_FILE | tr '\r\n' ' ')
perm_training_files_peas=$(awk -F "\"*,\"*" '{print $9;}' $PARAM_FILE | tr '\r\n' ' ')
perm_ground_truth_peas=$(awk -F "\"*,\"*" '{print $10;}' $PARAM_FILE | tr '\r\n' ' ')
wig_peas=$(awk -F "\"*,\"*" '{print $11;}' $PARAM_FILE | tr '\r\n' ' ')

# Parameters for all algorithms.
grid_size = 100
epochs = 100
epochs_cagt = 1000
cell_type_count=${#cell_types[@]}

# Run for each chromosome
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    # Repeat each run 5 times for robustness.
    for repeat in 1 2 3 4 5
    do
        # Make training and testing split directories.
        if [[ ! -e ${training_files[$i]}/shifted_split_training_$repeat ]]; then
            mkdir -p ${training_files[$i]}/shifted_split_training_$repeat
        fi
        if [[ ! -e ${training_files[$i]}/shifted_split_training_$repeat ]]; then
            mkdir -p ${training_files[$i]}/shifted_split_training_$repeat
        fi
        if [[ ! -e ${wig}_split_training_$repeat ]]; then
            mkdir -p ${wig}_split_training_$repeat
        fi
        if [[ ! -e ${wig}_split_testing_$repeat ]]; then
            mkdir -p ${wig}_split_testing_$repeat
        fi
        # if [[ ! -e ${perm_training_files[$i]}/shifted_split_training_$repeat ]]; then
            # mkdir -p ${perm_training_files[$i]}/shifted_split_training_$repeat
        # fi
        # if [[ ! -e ${perm_training_files[$i]}/shifted_split_testing_$repeat ]]; then
            # mkdir -p ${perm_training_files[$i]}/shifted_split_testing_$repeat
        # fi
        # if [[ ! -e ${training_files_peas[$i]}/shifted_split_training_$repeat ]]; then
            # mkdir -p ${training_files_peas[$i]}/shifted_split_training_$repeat
        # fi
        # if [[ ! -e ${training_files_peas[$i]}/shifted_split_testing_$repeat ]]; then
            # mkdir -p ${training_files_peas[$i]}/shifted_split_testing_$repeat
        # fi
        # if [[ ! -e ${wig_peas}_split_training_$repeat ]]; then
            # mkdir -p ${wig_peas}_split_training_$repeat
        # fi
        # if [[ ! -e ${wig_peas}_split_testing_$repeat ]]; then
            # mkdir -p ${wig_peas}_split_testing_$repeat
        # fi
        # if [[ ! -e ${perm_training_files_peas[$i]}/shifted_split_training_$repeat ]]; then
            # mkdir -p ${perm_training_files_peas[$i]}/shifted_split_training_$repeat
        # fi
        # if [[ ! -e ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat ]]; then
            # mkdir -p ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat
        # fi
        
        # Run for each cell type.
        for (( i=1; i<${cell_type_count}+1; i++ ));
        do
            
            # Randomly choose half of the regions for training.
            python random_split.py ${training_files[$i]}/shifted/${chrom}.pkl ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${training_files[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
            
            python random_split.py ${perm_training_files[$i]}/shifted/${chrom}.pkl ${perm_training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${perm_training_files[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
            
            python random_split.py ${training_files_peas[$i]}/shifted/${chrom}.pkl ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${training_files_peas[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
            
            python random_split.py ${perm_training_files_peas[$i]}/shifted/${chrom}.pkl ${perm_training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl ${perm_training_files_peas[$i]}/shifted_split_testing_$repeat/${chrom}.pkl
            
            python random_split_wig.py ${wig}/$chrom.wig ${wig}_split_testing_$repeat/$chrom.wig ${wig}_split_training_$repeat/$chrom.wig
            
            python random_split_wig.py ${wig_peas}/$chrom.wig ${wig_peas}_split_testing_$repeat/$chrom.wig ${wig_peas}_split_training_$repeat/$chrom.wig
        
            # Submit all SOM-VN and SOM jobs.
            for alpha_0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
            do
                for sigma_0 in 1 2 3 4 5 6 7 8 9
                do
                    # Create a directory.
                    if [[ ! -e ${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_somvn_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_somvn_training/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_somvn_signalperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_somvn_signalperm/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_somvn_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_somvn_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_somvn_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_somvn_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_somvn_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_somvn_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_som/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_som/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    if [[ ! -e ${cell_types[$i]}_som_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                        mkdir -p ${cell_types[$i]}_som_training/$chrom/$repeat/$alpha_0/$sigma_0
                    fi
                    
                    # Normal SOM-VN
                    qsub -A $PROJECT -v a=${training_files[$i]}/shifted/${chrom}.pkl -v b=$BIN_SIZE -v c=$chrom -v d=${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0/ -v g=$grid_size -v h=${chromhmm_anno[$i]} -v i=$epochs -v l=$alpha_0 -v n=$sigma_0 -v r=$REGION_SIZE -v s=$SCRIPTS -v t=$CCCUTOFF -v u=${training_files[$i]}/percentile_cutoffs/${chrom}.txt -v z=False learn_shapes_for_chrom_som_vn.sh

                    qsub -A $PROJECT -v a=${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl -v b=$BIN_SIZE -v c=$chrom -v d=${cell_types[$i]}_somvn/$chrom/$repeat/$alpha_0/$sigma_0/ -v g=$grid_size -v h=${chromhmm_anno[$i]} -v i=$epochs -v l=$alpha_0 -v n=$sigma_0 -v r=$REGION_SIZE -v s=$SCRIPTS -v t=$CCCUTOFF -v u=${training_files[$i]}/percentile_cutoffs/${chrom}.txt -v z=False learn_shapes_for_chrom_som_vn.sh                                         

                    # Permuted WIG SOM-VN
                    qsub -A $PROJECT -v a=${perm_training_files[$i]}/shifted/${chrom}.pkl -v b=$BIN_SIZE -v c=$chrom -v d=${cell_types[$i]}_somvn_signalperm/$chrom/$repeat/$alpha_0/$sigma_0/ -v g=$grid_size -v h=${chromhmm_anno[$i]} -v i=$epochs -v l=$alpha_0 -v n=$sigma_0 -v r=$REGION_SIZE -v s=$SCRIPTS -v t=$CCCUTOFF -v u=${perm_training_files[$i]}/percentile_cutoffs/${chrom}.txt -v z=False learn_shapes_for_chrom_som_vn.sh
                    
                    qsub -A $PROJECT -v a=${perm_training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl -v b=$BIN_SIZE -v c=$chrom -v d=${cell_types[$i]}_somvn_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0/ -v g=$grid_size -v h=${chromhmm_anno[$i]} -v i=$epochs -v l=$alpha_0 -v n=$sigma_0 -v r=$REGION_SIZE -v s=$SCRIPTS -v t=$CCCUTOFF -v u=${perm_training_files[$i]}/percentile_cutoffs/${chrom}.txt -v z=False learn_shapes_for_chrom_som_vn.sh
                    
                    # Permuted annotation SOM-VN
                    qsub -A $PROJECT -v a=${training_files[$i]}/shifted/${chrom}.pkl -v b=$BIN_SIZE -v c=$chrom -v d=${cell_types[$i]}_somvn_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0/ -v g=$grid_size -v h=${perm_chromhmm_anno[$i]} -v i=$epochs -v l=$alpha_0 -v n=$sigma_0 -v r=$REGION_SIZE -v s=$SCRIPTS -v t=$CCCUTOFF -v u=${training_files[$i]}/percentile_cutoffs/${chrom}.txt -v z=False learn_shapes_for_chrom_som_vn.sh
                    
                    qsub -A $PROJECT -v a=${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl -v b=$BIN_SIZE -v c=$chrom -v d=${cell_types[$i]}_somvn_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0/ -v g=$grid_size -v h=${perm_chromhmm_anno[$i]} -v i=$epochs -v l=$alpha_0 -v n=$sigma_0 -v r=$REGION_SIZE -v s=$SCRIPTS -v t=$CCCUTOFF -v u=${training_files[$i]}/percentile_cutoffs/${chrom}.txt -v z=False learn_shapes_for_chrom_som_vn.sh
                    
                    # Normal SOM
                    qsub -A $PROJECT learn_shapes_for_chrom_som.sh -a ${training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_som/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${chromhmm_anno[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z False
                    
                    qsub -A $PROJECT learn_shapes_for_chrom_som.sh -a ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_som_training/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${chromhmm_anno[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z False
                    
                    # Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} = "GM12878"]
                    then
                    
                        if [[ ! -e ${cell_types[$i]}_somvn_peas/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_somvn_peas/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_somvn_peas_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_somvn_peas_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_somvn_peas_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_somvn_peas_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_somvn_peas_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_somvn_peas_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_somvn_peas_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_somvn_peas_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_som_peas/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_som_peas/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        if [[ ! -e ${cell_types[$i]}_som_peas_training/$chrom/$repeat/$alpha_0/$sigma_0 ]]; then
                            mkdir -p ${cell_types[$i]}_som_peas_training/$chrom/$repeat/$alpha_0/$sigma_0
                        fi
                        # # Normal SOM-VN with PEAS ground truth
                        # qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_somvn_peas/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt -z True
                        
                        # qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_somvn_peas_training/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt -z True
                        
                        # # Permuted WIG SOM-VN with PEAS ground truth
                        # qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${perm_training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_somvn_peas_signalperm/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${perm_training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt -z True

                        # qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${perm_training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_somvn_peas_signalperm_training/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${perm_training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt -z True
                        
                        # # Permuted annotation SOM-VN with PEAS ground truth
                        # qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_somvn_peas_chromhmmperm/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${perm_ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt -z True
                        
                        # qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_somvn_peas_chromhmmperm_training/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${perm_ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.txt -z True
                        
                        # # Normal SOM with PEAS ground truth
                        # qsub -A $PROJECT learn_shapes_for_chrom_som.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_som_peas/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z True
                        
                        # qsub -A $PROJECT learn_shapes_for_chrom_som.sh -a ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_som_peas_training/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z True
                    # fi
                done
            done
            
            # Submit all CAGT jobs.
            for k in 20 40 80 160
            do
                for max_dist in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
                do
                    
                    if [[ ! -e ${cell_types[$i]}_cagt/$chrom/$repeat/$k/$max_dist ]]; then
                            mkdir -p ${cell_types[$i]}_cagt/$chrom/$repeat/$k/$max_dist
                    fi
                    if [[ ! -e ${cell_types[$i]}_cagt_training/$chrom/$repeat/$k/$max_dist ]]; then
                            mkdir -p ${cell_types[$i]}_cagt_training/$chrom/$repeat/$k/$max_dist
                    fi
                    
                    # CAGT
                    qsub -A $PROJECT learn_shapes_for_chrom_cagt.sh -a ${training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_cagt/$chrom/$repeat/$k/$max_dist/ -h ${chromhmm_anno[$i]} -i $epochs_cagt -k $k -m $max_dist -p $CAGT_PATH -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z False
                    
                    qsub -A $PROJECT learn_shapes_for_chrom_cagt.sh -a ${training_files[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_cagt_training/$chrom/$repeat/$k/$max_dist/ -h ${chromhmm_anno[$i]} -i $epochs_cagt -k $k -m $max_dist -p $CAGT_PATH -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z False
                    
                    # # Do the same for PEAS ground truth for GM12878.
                    # if [ ${cell_types[$i]} = "GM12878"]
                    # then
                    
                        if [[ ! -e ${cell_types[$i]}_cagt_peas/$chrom/$repeat/$k/$max_dist ]]; then
                            mkdir -p ${cell_types[$i]}_cagt_peas/$chrom/$repeat/$k/$max_dist
                        fi
                        if [[ ! -e ${cell_types[$i]}_cagt_peas_training/$chrom/$repeat/$k/$max_dist ]]; then
                                mkdir -p ${cell_types[$i]}_cagt_peas_training/$chrom/$repeat/$k/$max_dist
                        fi
                    
                        # # CAGT with PEAS ground truth
                        # qsub -A $PROJECT learn_shapes_for_chrom_cagt.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_cagt_peas/$chrom/$repeat/$k/$max_dist/ -h ${training_files_peas[$i]} -i $epochs_cagt -k $k -m $max_dist -p $CAGT_PATH -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z True
                        
                        # qsub -A $PROJECT learn_shapes_for_chrom_cagt.sh -a ${training_files_peas[$i]}/shifted_split_training_$repeat/${chrom}.pkl -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_cagt_peas_training/$chrom/$repeat/$k/$max_dist/ -h ${training_files_peas[$i]} -i $epochs_cagt -k $k -m $max_dist -p $CAGT_PATH -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -z True
                    # fi
                done
            done
            
            # Submit all intensity-only jobs.
            
            if [[ ! -e ${cell_types[$i]}_signal/$chrom/$repeat ]]; then
                    mkdir -p ${cell_types[$i]}_signal/$chrom/$repeat
            fi
            if [[ ! -e ${cell_types[$i]}_signal_training/$chrom/$repeat ]]; then
                    mkdir -p ${cell_types[$i]}_signal_training/$chrom/$repeat
            fi
            
            qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_signal/$chrom/$repeat/ -h ${chromhmm_anno[$i]} -s $SCRIPTS -w $wig/$chrom.wig -z False
            
            qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_signal_training/$chrom/$repeat/ -h ${chromhmm_anno[$i]} -s $SCRIPTS -w ${wig}_split_training_$repeat/$chrom.wig -z False
            
            # Submit all intensity-only jobs with PEAS ground truth.
            # if [ ${cell_types[$i]} = "GM12878"]
            # then
            
                # if [[ ! -e ${cell_types[$i]}_peas_signal/$chrom/$repeat ]]; then
                    # mkdir -p ${cell_types[$i]}_peas_signal/$chrom/$repeat
                # fi
                # if [[ ! -e ${cell_types[$i]}_peas_signal_training/$chrom/$repeat ]]; then
                        # mkdir -p ${cell_types[$i]}_peas_signal_training/$chrom/$repeat
                # fi
            
                # qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_peas_signal/$chrom/$repeat/ -h ${ground_truth_peas[$i]} -s $SCRIPTS -w $wig_peas/$chrom.wig -z True
                
                # qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -b $BIN_SIZE -c $chrom -d ${cell_types[$i]}_peas_signal_training/$chrom/$repeat/ -h ${ground_truth_peas[$i]} -s $SCRIPTS -w ${wig_peas}_split_training_$repeat/$chrom.wig -z True
            # fi
            
        done
    done
done
