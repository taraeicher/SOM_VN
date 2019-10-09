#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   

#Variables
    USAGE="This script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using CAGT, then merged to correct for signal shift. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-i> The bin size used to generate the WIG file (default: 50 bp)\n
    <-r> The size of the input regions (default: 4000)\n
    <-t> The cutoff to use for cross-correlation significance.\n
    <-a> Directory containing training regions\n
    <-u> Percentile cutoff file\n
    <-c> The chromosome name\n
    <-p> The path to the CAGT file\n
    <-s> The directory containing the scripts"
    
    echo -e $USAGE
    REGION_SIZE=4000
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=50
    TRAINING=""
    CCCUTOFF=0.75
    CAGT_PATH="/fs/project/PAS0272/Tara/DNase_SOM/Brain/cagt/matlab/src"
    SCRIPTS=""
    while getopts h:d:c:i:r:t:a:u:p: option; do
        case "${option}" in
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            a) TRAINING=$(realpath $OPTARG);;
            u) CUTOFFS=$(realpath $OPTARG);;
            c) CHROM=$(realpath $OPTARG);;
            t) CCCUTOFF=$OPTARG;;
            p) CAGT_PATH=$(realpath $OPTARG);;
            s) SCRIPTS=$(realpath $OPTARG);;
        esac
    done
	
    #Move to the directory containing the scripts.
    cd $SCRIPTS
    
    #Create all needed directories.
    MATLAB_MATRIX="$BASE_PATH/cagt_matlab_mat"
    if [[ ! -e $MATLAB_MATRIX ]]; then
        mkdir $MATLAB_MATRIX
    fi
    CAGT_OUT="$BASE_PATH/cagt_output"
    if [[ ! -e $CAGT_OUT ]]; then
        mkdir $CAGT_OUT
    fi
    CAGT_OUT_CSV="$BASE_PATH/cagt_output_csv"
    if [[ ! -e $CAGT_OUT_CSV ]]; then
        mkdir $CAGT_OUT_CSV
    fi
    CAGT_SHIFTED="$BASE_PATH/cagt_output_shifted"
    if [[ ! -e $CAGT_SHIFTED ]]; then
        mkdir $CAGT_SHIFTED
    fi
    ANNOTATED="$BASE_PATH/cagt_anno_beds"
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
    INTERSECTS="$BASE_PATH/cagt_intersects"
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    INTERSECTS_SORTED="$BASE_PATH/cagt_intersects_sorted"
    if [[ ! -e $INTERSECTS_SORTED ]]; then
        mkdir $INTERSECTS_SORTED
    fi
    CHROMHMM_DISTRIB="$BASE_PATH/cagt_chromhmm_distrib"
    if [[ ! -e $CHROMHMM_DISTRIB ]]; then
        mkdir $CHROMHMM_DISTRIB
    fi

    #Extract the signal and run CAGT.
    python extract_signal.py $TRAINING $TRAINING_CSV/$CHROM.csv
    module load matlab
    matlab -nodisplay "run_cagt_func $TRAINING_CSV/$CHROM.csv $MATLAB_MATRIX/$CHROM.mat $CAGT_OUT_CSV/$CHROM.csv $CAGT_PATH"
    matlab -nodisplay -nodesktop -r "run_cagt('$TRAINING_CSV/$CHROM.csv','$MATLAB_MATRIX/$CHROM.mat','$CAGT_OUT_CSV/$CHROM.csv','$CAGT_PATH')"
    python convert_to_pickle.py $CAGT_OUT_CSV $CAGT_OUT
    echo -e "SOM model is ready for chrom $CHROM.\n"
    
    #Merge shifted regions.
    python merge_shifted.py $CAGT_OUT/$CHROM.pkl $CAGT_SHIFTED/$CHROM.pkl $CCCUTOFF
    echo -e "Merging complete for chrom $CHROM.\n"
    
    #Annotate regions with shape.
    python make_shape_bed.py $TRAINING/$CHROM.pkl $CAGT_SHIFTED/$CHROM.pkl $ANNOTATED/$CHROM.bed
    echo -e "Initial annotations complete for chrom $CHROM.\n"
    
    #Intersect regions with ChromHMM.
    bedtools intersect -wao -a $ANNOTATED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $SOM_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl
        
    #Exit
	wait
	echo Done!
	exit 0