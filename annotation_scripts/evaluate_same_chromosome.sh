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

promoter_cutoff=0.05
enhancer_cutoff=0.05
repressor_cutoff=0.9
weak_cutoff=0.9

cell_type_count=${#cell_types[@]}
for chrom in 1 8 15 22
do
    for (( i=1; i<${cell_type_count}+1; i++ ));
    do
        for repeat in 1 2 3 4 5
        do
            for alpha_0 in 0.2 0.4 0.6 0.8
            do
                for sigma_0 in 2 4 6 8
                do
                
                    som_vn_basepath=$BASE_PATH/som_vn_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                    som_vn_perm_basepath=$BASE_PATH/som_vn_signalperm_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                    som_vn_chromhmmperm_basepath=$BASE_PATH/som_vn_chromhmmperm_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                    som_basepath=$BASE_PATH/som_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                    
                    if [[ ! -e $som_vn_basepath ]]; then
                        mkdir -p $som_vn_basepath
                    fi
                    if [[ ! -e $som_vn_perm_basepath ]]; then
                        mkdir -p $som_vn_perm_basepath
                    fi
                    if [[ ! -e $som_vn_chromhmmperm_basepath ]]; then
                        mkdir -p $som_vn_chromhmmperm_basepath
                    fi
                    if [[ ! -e $som_basepath ]]; then
                        mkdir -p $som_basepath
                    fi
                    
                    #SOM-VN Normal
                    qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_vn_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_vn/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=False
                   
                    #SOM-VN Permuted Signal
                    qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_vn_perm_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_vn_signalperm/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=False
                   
                    #SOM-VN Permuted ChromHMM
                    qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_vn_chromhmmperm_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_vn_chromhmmperm/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=False
                   
                    #SOM
                    qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=False
                   
                    if [ ${cell_types[$i]} = "GM12878"];
                    then
                        som_vn_peas_basepath=$BASE_PATH/som_vn_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                        som_vn_peas_perm_basepath=$BASE_PATH/som_vn_signalperm_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                        som_vn_peas_chromhmmperm_basepath=$BASE_PATH/som_vn_chromhmmperm_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                        som_peas_basepath=$BASE_PATH/som_peas_intrachrom/${cell_types[$i]}/$chrom/$repeat/$alpha_0/$sigma_0/
                        
                        if [[ ! -e $som_vn_peas_basepath ]]; then
                            mkdir -p $som_vn_peas_basepath
                        fi
                        if [[ ! -e $som_vn_peas_perm_basepath ]]; then
                            mkdir -p $som_vn_peas_perm_basepath
                        fi
                        if [[ ! -e $som_vn_peas_chromhmmperm_basepath ]]; then
                            mkdir -p $som_vn_peas_chromhmmperm_basepath
                        fi
                        if [[ ! -e $som_peas_basepath ]]; then
                            mkdir -p $som_peas_basepath
                        fi
                    
                        #SOM-VN PEAS
                        qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_vn_peas_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno_peas[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_vn_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=True

                        #SOM-VN Permuted PEAS                       
                        qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_vn_peas_perm_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno_peas[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_vn_signalperm_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=True
                   
                        #SOM-VN ChromHMM Permuted PEAS
                        qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_vn_peas_chromhmmperm_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno_peas[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_vn_chromhmmperm_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=True
                   
                        #SOM PEAS
                        qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $som_peas_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno_peas[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/som_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$alpha_0/$sigma_0/) REGIONS=$(realpath ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=True
                    fi
                done
            done
            # Submit all CAGT jobs.
            for k in 20 40 80 160
            do
                for max_dist in 0.3 0.5 0.7 0.9
                do
                    cagt_basepath=$BASE_PATH/cagt/${cell_types[$i]}_split_training/$chrom/$repeat/$k/$max_dist/
                    
                    if [[ ! -e $cagt_basepath ]]; then
                        mkdir -p $cagt_basepath
                    fi
                    
                    # CAGT                   
                    qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $cagt_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_anno[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/cagt/${cell_types[$i]}_split_training/$chrom/$repeat/$k/$max_dist/) REGIONS=$(realpath ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=False

                    # Do the same for PEAS ground truth for GM12878.
                    if [ ${cell_types[$i]} = "GM12878"];
                    then
                        cagt_peas_basepath=$BASE_PATH/cagt_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$k/$max_dist/
                    
                        if [[ ! -e $cagt_peas_basepath ]]; then
                            mkdir -p $cagt_peas_basepath
                        fi
                        
                        # CAGT with PEAS ground truth
                        qsub -A $PROJECT annotate_and_get_pr.sh PEAKS=$(realpath ${peaks[$i]}) CHROM=$chrom BASE_PATH=$(realpath $cagt_peas_basepath) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff CHROMHMM=$(realpath ${chromhmm_peas_anno[$i]}) SCRIPTS=$(realpath $SCRIPTS) SHAPES=$(realpath $BASE_PATH/cagt_peas/${cell_types[$i]}_split_training/$chrom/$repeat/$k/$max_dist/) REGIONS=$(realpath ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl) IS_PEAS=True
                    fi
                done
            done
            
            signal_basepath=$BASE_PATH/signal_intrachrom/${cell_types[$i]}/$chrom/
                    
            if [[ ! -e $signal_basepath ]]; then
                mkdir -p $signal_basepath
            fi
            
            #Signal only
            qsub -A $PROJECT annotate_and_get_pr_signal.sh PEAKS=${peaks[$i]} REGIONS=$(realpath ${training_files[$i]}/shifted_split_testing_$repeat/$chrom.pkl) SIGNALS=$(realpath $BASE_PATH/signal/${cell_types[$i]}_split_training/$chrom/) BASE_PATH=$(realpath $signal_basepath) CHROM=$chrom CHROMHMM=$(realpath ${chromhmm_anno[$i]}) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff SCRIPTS=$(realpath $SCRIPTS) IS_PEAS=False
            
            if [ ${cell_types[$i]} = "GM12878"]
            then
                signal_peas_basepath=$BASE_PATH/signal_peas_intrachrom/${cell_types[$i]}/$chrom/
                    
                if [[ ! -e $signal_peas_basepath ]]; then
                    mkdir -p $signal_peas_basepath
                fi
                
                #Signal only PEAS
                qsub -A $PROJECT annotate_and_get_pr_signal.sh PEAKS=${peaks[$i]} REGIONS=$(realpath ${training_files_peas[$i]}/shifted_split_testing_$repeat/$chrom.pkl) SIGNALS=$(realpath $BASE_PATH/signal_peas/${cell_types[$i]}_split_training/$chrom/) BASE_PATH=$(realpath $signal_peas_basepath) CHROM=$chrom CHROMHMM=$(realpath ${chromhmm_anno_peas[$i]}) PROMOTER_CUTOFF=$promoter_cutoff ENHANCER_CUTOFF=$enhancer_cutoff REPRESSOR_CUTOFF=$repressor_cutoff WEAK_CUTOFF=$weak_cutoff SCRIPTS=$(realpath $SCRIPTS) IS_PEAS=True
            fi
        done
    done
done
