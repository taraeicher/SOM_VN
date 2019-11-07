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
cell_types=( $(cut -d ',' -f1 $PARAM_FILE ) )
training_files=( $(cut -d ',' -f2 $PARAM_FILE ) )
perm_training_files=( $(cut -d ',' -f3 $PARAM_FILE ) )
chromhmm_anno=( $(cut -d ',' -f4 $PARAM_FILE ) )
perm_chromhmm_anno=( $(cut -d ',' -f5 $PARAM_FILE ) )
wig=( $(cut -d ',' -f6 $PARAM_FILE ) )
training_files_peas=( $(cut -d ',' -f7 $PARAM_FILE ) )
ground_truth_peas=( $(cut -d ',' -f8 $PARAM_FILE ) )
perm_training_files_peas=( $(cut -d ',' -f9 $PARAM_FILE ) )
perm_ground_truth_peas=( $(cut -d ',' -f10 $PARAM_FILE ) )
wig_peas=( $(cut -d ',' -f11 $PARAM_FILE ) )

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
        # Randomly choose half of the regions for training.
        
        # Run for each cell type.
        for (( i=1; i<${cell_type_count}+1; i++ ));
        do
            # Submit all SOM-VN and SOM jobs.
            for alpha_0 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
            do
                for sigma_0 in 1 2 3 4 5 6 7 8 9
                do
                    # Normal SOM-VN
                    qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_vn/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${chromhmm_anno[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files[$i]}/percentile_cutoffs/${chrom}.pkl     

                    # Permuted WIG SOM-VN
                    qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${perm_training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_vn_signalperm/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${chromhmm_anno[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${perm_training_files[$i]}/percentile_cutoffs/${chrom}.pkl 
                    
                    # Permuted annotation SOM-VN
                    qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_vn_chromhmmmperm/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${perm_chromhmm_anno[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files[$i]}/percentile_cutoffs/${chrom}.pkl 
                    
                    # Normal SOM
                    qsub -A $PROJECT learn_shapes_for_chrom_som.sh -a ${training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${chromhmm_anno[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF
                    
                    # Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} = "GM12878"]
                    then
                        # Normal SOM-VN with PEAS ground truth
                        qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_vn_peas/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.pkl 
                        
                        # Permuted WIG SOM-VN with PEAS ground truth
                        qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${perm_training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_vn_signalperm_peas/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${perm_training_files_peas[$i]}/percentile_cutoffs/${chrom}.pkl 
                        
                        # Permuted annotation SOM-VN with PEAS ground truth
                        qsub -A $PROJECT learn_shapes_for_chrom_som_vn.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_vn_chromhmmmperm_peas/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${perm_ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF -u ${training_files_peas[$i]}/percentile_cutoffs/${chrom}.pkl 
                        
                        # Normal SOM with PEAS ground truth
                        qsub -A $PROJECT learn_shapes_for_chrom_som.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/som_peas/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/ -g $grid_size -h ${ground_truth_peas[$i]} -i $epochs -l $alpha_0 -n $sigma_0 -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF
                    fi
                done
            done
            
            # Submit all CAGT jobs.
            for k in 20 40 80 160
            do
                for max_dist in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
                do
                    # CAGT
                    qsub -A $PROJECT learn_shapes_for_chrom_cagt.sh -a ${training_files[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/cagt/${cell_types[$i]}/$chrom/$repeat/$k/$max_dist/ -h ${chromhmm_anno[$i]} -i $epochs_cagt -k $k -m $max_dist -p $CAGT_PATH -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF
                    
                    # Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} = "GM12878"]
                    then
                        # CAGT with PEAS ground truth
                        qsub -A $PROJECT learn_shapes_for_chrom_cagt.sh -a ${training_files_peas[$i]}/shifted/${chrom}.pkl -b $BIN_SIZE -c $chrom -d $BASE_PATH/cagt_peas/${cell_types[$i]}/$chrom/$repeat/$k/$max_dist/ -h ${training_files_peas[$i]} -i $epochs_cagt -k $k -m $max_dist -p $CAGT_PATH -r $REGION_SIZE -s $SCRIPTS -t $CCCUTOFF
                    fi
                done
            done
            
            # Submit all intensity-only jobs.
            qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -b $BIN_SIZE -c $chrom -d $BASE_PATH/signal/${cell_types[$i]}/$chrom/$repeat/ -h ${chromhmm_anno[$i]} -s $SCRIPTS -w $wig/$chrom.wig
            
            # Submit all intensity-only jobs with PEAS ground truth.
            if [ ${cell_types[$i]} = "GM12878"]
            then
                qsub -A $PROJECT learn_shapes_for_chrom_signal.sh -b $BIN_SIZE -c $chrom -d $BASE_PATH/signal_peas/${cell_types[$i]}/$chrom/$repeat/ -h ${ground_truth_peas[$i]} -s $SCRIPTS -w $wig_peas/$chrom.wig
            fi
            
        done
    done
done
