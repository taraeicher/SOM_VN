Skip to content
Search or jump to…

Pull requests
Issues
Marketplace
#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   

#Variables
    USAGE="\n\nThis script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using CAGT, then merged to correct for signal shift. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n\n
    <-a> Directory containing training regions\n
    <-b> The bin size used to generate the WIG file (default: 10 bp)\n
    <-c> The chromosome name\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-h> The ChromHMM file used for intersecting.\n
    <-i> The maximum number of iterations for CAGT\n
    <-k> The number of clusters to learn prior to agglomerative clustering\n
    <-m> The maximum distance for merging to occur in the agglomerative clustering step of CAGT (default: 0.8)\n
    <-p> The path to the CAGT file\n
    <-r> The size of the input regions (default: 1000)\n
    <-s> The directory containing the scripts\n
    <-t> cross-correlation cutoff\n
    <-z> Whether or not this we are using PEAS annotations\n\n"
    
    echo -e $USAGE
    REGION_SIZE=1000
    BASE_PATH=""
    BAM=""
    CHROMHMM=""
    BIN_SIZE=10
    TRAINING=""
    CCCUTOFF=0.75
    CAGT_PATH="/fs/project/PAS0272/Tara/DNase_SOM/Brain/cagt/matlab/src"
    SCRIPTS=""
    MERGE_DIST=0.8
    ITERATIONS=1000
    K=40
    IS_PEAS=False
    while getopts a:b:c:d:h:i:k:m:p:r:s:t:z: option; do
        case "${option}" in
            a) TRAINING=$(realpath $OPTARG);;
            b) BIN_SIZE=$OPTARG;;
            c) CHROM=$OPTARG;;
            d) BASE_PATH=$(realpath $OPTARG);;
            h) CHROMHMM=$(realpath $OPTARG);;
            i) ITERATIONS=$OPTARG;;
            k) K=$OPTARG;;
            m) MERGE_DIST=$OPTARG;;
            p) CAGT_PATH=$(realpath $OPTARG);;
            r) REGION_SIZE=$OPTARG;;
            s) SCRIPTS=$(realpath $OPTARG);;
            t) CCCUTOFF=$OPTARG;;
            z) IS_PEAS=$OPTARG;;
        esac
    done
	
    #Move to the directory containing the scripts.
    cd $SCRIPTS
    
    #Create all needed directories.
    TRAINING_CSV="$BASE_PATH/training_csv"
    if [[ ! -e $TRAINING_CSV ]]; then
        mkdir $TRAINING_CSV
    fi
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
    #python extract_signal.py $TRAINING $TRAINING_CSV/$CHROM.csv
    module load matlab
    matlab -nodisplay -nodesktop -r "run_cagt('$TRAINING/$CHROM,'$MATLAB_MATRIX/$CHROM.mat','$CAGT_OUT_CSV/$CHROM.csv','$CAGT_PATH', '$MERGE_DIST', '$ITERATIONS', '$K')"
	matlab -nodisplay -nodesktop -r "run_cagt('/data/eichertd/som_vn_data/training_b_cell_low_shifted/chrom1','/data/eichertd/som_vn_data/matlab_matrix_b_cell_low/','1','/data/eichertd/som_vn_data/cagt_out_b_cell_low/chrom1.csv','/home/eichertd/som_vn_code/SOM_VN/cagt/trunk/matlab/src/', '0.8', '1000', '40')"
    python convert_to_pickle.py $CAGT_OUT_CSV/$CHROM.csv $CAGT_OUT/$CHROM.pkl
    echo -e "CAGT model is ready for chrom $CHROM.\n"
    
    #Merge shifted regions.
    python merge_shifted.py $CAGT_OUT/$CHROM.pkl $CAGT_SHIFTED/$CHROM.pkl $CCCUTOFF
    echo -e "Merging complete for chrom $CHROM.\n"
    
    #Annotate regions with shape.
    python make_shape_bed.py $TRAINING $CAGT_SHIFTED/$CHROM.pkl $ANNOTATED/$CHROM.bed
    echo -e "Initial annotations complete for chrom $CHROM.\n"
    
    #Intersect regions with ChromHMM.
    bedtools intersect -wao -a $ANNOTATED/$CHROM.bed -b $CHROMHMM > $INTERSECTS/$CHROM.bed
    bedtools sort -i $INTERSECTS/$CHROM.bed > $INTERSECTS_SORTED/$CHROM.bed
    python find_chromhmm_distrib.py $INTERSECTS_SORTED/$CHROM.bed $CAGT_SHIFTED/$CHROM.pkl $CHROMHMM_DISTRIB/$CHROM.pkl $IS_PEAS
        
    #Exit
	wait
	echo Done!
	exit 0
