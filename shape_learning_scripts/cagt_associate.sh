#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   
# """
# Dependencies:
# 1. Following python modules for gosr: argparse, logging, sys, itertools, numpy, pysam
# 2. Tensorflow, which can be installed here: https://www.tensorflow.org/install/install_linux
# 3. Samtools: https://sourceforge.net/projects/samtools/files/
# 4. Wig-split: https://github.com/taoliu/taolib/blob/master/Scripts/wig_split.py
# 5. Taolib directory: https://github.com/taoliu/taolib
# 6. Gap statistic code: https://github.com/minddrummer/gap
# 7. Added to your path: gosr/bin, taolib, wig_split, gap
# 8. Added to your python path: pythonmodules/gosr/lib
# 9. Make sure you are not using python3, which causes an error with the print statements in gosr and in wig_split
# 10. Make sure your version of python is compiled to use UCS2

#Variables
    CELL_LINE=""
    REGION_SIZE=4000
    BASE_FILENAME=""
    BAM=""
    CHROMHMM=""
    CHROMS_NUM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    WIG_SPLIT_PATH=""
    BIN_SIZE=50
    SHAPES_COMPREHENSIVE=""
    while getopts n:d:c:i:r:s: option; do
        case "${option}" in
            n) CELL_LINE=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            c) CHROMHMM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            s) SHAPES_COMPREHENSIVE=$(realpath $OPTARG);;
        esac
    done
    BASE_PATH=$BASE_FILENAME/$CELL_LINE
    PYTHON_VERSION=2.7	
    
#Print message to user.
    echo "Annotating with the following settings:"
    echo "Name is $CELL_LINE"
    echo "Base directory is $BASE_FILENAME"
    echo "ChromHMM file is $CHROMHMM"
    echo "Bin size in WIG file is $BIN_SIZE"
    echo "Region size for annotation is $REGION_SIZE"
    echo "Shapes will be saved in $SHAPES_COMPREHENSIVE"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"
	
#Create all needed directories.
    CAGT_OUT="$BASE_PATH/cagt/matlab/src"
    if [[ ! -e $CAGT_OUT ]]; then
        mkdir $CAGT_OUT
    fi
    TRAINING_ANNOTATION_FILES="$BASE_PATH/training_anno_allchrom"
    if [[ ! -e $TRAINING_ANNOTATION_FILES ]]; then
        mkdir $TRAINING_ANNOTATION_FILES
    fi
    WIG="$BASE_PATH/wig_chroms"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    SHAPE_ANNOTATED="$BASE_PATH/anno_beds_cagt"
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED="$BASE_PATH/anno_beds_sorted_cagt"
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL="$BASE_PATH/anno_beds_final_cagt"
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects_cagt"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    if [[ ! -e ${SHAPES_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}_${CELL_LINE}
    fi
    SHAPES="$BASE_PATH/shapes_cagt_merged"
    SHAPES_LOG="$BASE_PATH/shapes_log_cagt"
	
# Method for running the pipeline for a chromocagte.
	run_pipeline() {
		local c=$1
		
		# #Annotate regions with shape.
		# python make_shape_bed.py $TRAINING_ANNOTATION_FILES/chrom${c}window3 $CAGT_OUT/ctcf_centroids_${c}.csv $SHAPE_ANNOTATED/anno$c 0
        # echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
        
        # bedtools sort -i  $SHAPE_ANNOTATED/anno${c} > $SHAPE_ANNOTATED_SORTED/anno${c}
        # python consolidate.py $SHAPE_ANNOTATED_SORTED/anno${c} $SHAPE_ANNOTATED_FINAL/anno${c}
        # cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/anno${c}.bed
        # awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/clusters_anno${c}
        # cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/scores_anno${c}.bed
        # echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
        
        # #Save shapes to file.
        # bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${c}.bed -b $CHROMHMM/${CELL_LINE}_chromhmm_15_liftOver.bed > $CHROMHMM_INTERSECTS/anno${c}.bed			
        # bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
        python consolidate_chromHMM.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $CAGT_OUT/ctcf_centroids_${c}.csv $SHAPES_COMPREHENSIVE ${WIG}/${CELL_LINE}.chr${c}.wig ${c} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${c}window3 0
	}
    
    #Run the pipeline from split WIG files to final set of annotations.
    for f in $CHROMS_NUM;
        do 
            run_pipeline $f & 
        done
        
    #Merge shapes entries across chromocagtes.
    #python ../common_scripts/merge_significant.py $SHAPES_COMPREHENSIVE $SHAPES $SHAPES_LOG
        
    #Exit
	wait
	echo Done!
	exit 0