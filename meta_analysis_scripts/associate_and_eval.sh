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
    while getopts n:d:c:b:i:w:r:s: option; do
        case "${option}" in
            n) CELL_LINE=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            c) CHROMHMM=$(realpath $OPTARG);;
            b) BAM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            w) WIG_SPLIT_PATH=$(realpath $OPTARG);;
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
    echo "BAM file is $BAM"
    echo "Bin size in WIG file is $BIN_SIZE"
    echo "wig_split path is $WIG_SPLIT_PATH"
    echo "Region size for annotation is $REGION_SIZE"
    echo "Shapes will be saved in $SHAPES_COMPREHENSIVE"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"
	
#Create all needed directories
    SOM_OUT_FINAL="$BASE_PATH/som_output_final_allchrom"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG="$BASE_PATH/wig_chroms"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects_allchrom"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    TO_ANNOTATE="$BASE_PATH/annotation_files"
    if [[ ! -e $TO_ANNOTATE ]]; then
		mkdir $TO_ANNOTATE
	fi
    ANNOTATED="$BASE_PATH/annotated"
    if [[ ! -e $ANNOTATED ]]; then
		mkdir $ANNOTATED
	fi
    ANNOTATED_SORTED="$BASE_PATH/annotated_sorted"
    if [[ ! -e $ANNOTATED_SORTED ]]; then
		mkdir $ANNOTATED_SORTED
	fi
    ANNOTATED_CONSOLIDATED="$BASE_PATH/annotated_consolidated"
	if [[ ! -e $ANNOTATED_CONSOLIDATED ]]; then
		mkdir $ANNOTATED_CONSOLIDATED
	fi
    ANNO_MERGED="$BASE_PATH/annotated_merged"
	if [[ ! -e $ANNO_MERGED ]]; then
		mkdir $ANNO_MERGED
	fi
    SHAPES="$BASE_PATH/shapes"
    SHAPES_LOG="$BASE_PATH/shapes_log"
    WIGS="$BASE_PATH/wig_chroms_annotation"
    TRAINING_ANNOTATION_FILES="$BASE_PATH/training_anno_allchrom"
    PRECISION_RECALL="/fs/project/PAS0272/Tara/DNase_SOM/$CELL_LINE/pr_curve_$CELL_LINE"
    if [[ ! -e $PRECISION_RECALL ]]; then
		mkdir $PRECISION_RECALL
	fi

    anno_func() {
                local c=$1
                local cutoff=$2
                
                python ../annotation_scripts/make_annotated_bed.py $TO_ANNOTATE/chrom${c}window3 ${SHAPES}_${cutoff} $ANNOTATED/anno${c}_${cutoff} $WIGS/${CELL_LINE}.chr${c}.wig 0.0
                bedtools sort -i  $ANNOTATED/anno${c}_${cutoff} > $ANNOTATED_SORTED/anno${c}_${cutoff}
                bedtools sort -i  $ANNOTATED/anno${c}_${cutoff}clust > $ANNOTATED_SORTED/anno${c}_${cutoff}.clust
                python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED/anno${c}_${cutoff} $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff}
                cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff} > $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff}.bed
                cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff}.clust > $ANNOTATED_CONSOLIDATED/anno${c}clust_${cutoff}.bed
                awk '{ print $6}' $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff} > $ANNOTATED_CONSOLIDATED/clusters_anno${c}_${cutoff}
                cut -d$'\t' -f 7,8,9,10 $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff} > $ANNOTATED_CONSOLIDATED/scores_anno${c}_${cutoff}.bed
                
                # Intersect with ground truth.
                bedtools intersect -wao -a $ANNOTATED_CONSOLIDATED/anno${c}_${cutoff}.bed -b $CHROMHMM > $ANNO_MERGED/anno${c}_${cutoff}.bed
            }   
    
    #Run the pipeline from split WIG files to final set of annotations.
    cutoff_list="0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75"
    echo $cutoff_list > cutoffs.txt
    for cutoff in $cutoff_list;
        do
            #Create directory for output.
            if [[ ! -e ${SHAPES_COMPREHENSIVE}_${cutoff}_${CELL_LINE} ]]; then
                mkdir ${SHAPES_COMPREHENSIVE}_${cutoff}_${CELL_LINE}
            fi
            
            #Associate shapes with RE's for all chromosomes.
            # for c in $CHROMS_NUM;
                # do 
                    # python consolidate_chromHMM_cutoff.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $SOM_OUT_FINAL/chrom${c}som_centroid ${SHAPES_COMPREHENSIVE}_${cutoff} ${WIG}/${CELL_LINE}.chr${c}.wig ${c} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${c}window3 0  $cutoff
                # done
                
            #Merge shapes entries across chromosomes.
            #python ../common_scripts/merge_significant.py ${SHAPES_COMPREHENSIVE}_${cutoff} ${SHAPES}_${cutoff} ${SHAPES_LOG}_${cutoff}
            
            #Annotate
             # for c in $CHROMS_NUM;
                # do
                    # #Annotate
                    # anno_func $c $cutoff &
                # done
        done
        
    #Evaluate precision and recall.
    python plot_pr_curves.py $ANNO_MERGED/ $ANNOTATED_CONSOLIDATED/ $BASE_FILENAME/  $PRECISION_RECALL/ $WIG/$CELL_LINE.chr $CELL_LINE $CELL_LINE cutoffs.txt             
        
    #Exit
	wait
	echo Done!
	exit 0