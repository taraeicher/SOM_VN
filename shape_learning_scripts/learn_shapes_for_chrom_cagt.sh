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
    USAGE="This script is used for learning a set of representative shapes from the training regions output by create_regions.sh and annotating them with an RE. Shapes are learned for each chromosome using an SOM, then merged to correct for signal shift and clustered using k-means. Finally, the shapes are associated with RE by annotating the training set and associating the shapes with ChromHMM elements.\n
    <-n> The name of the cell line (e.g. Brain)\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-c> The ChromHMM file used for intersecting.\n
    <-i> The bin size used to generate the WIG file (default: 50 bp)\n
    <-r> The size of the input regions (default: 4000)\n
    <-s> The final set of shapes consolidated across all chromosomes.\n
    <-a> Directory containing training regions"
    
    echo -e $USAGE
    CELL_LINE=""
    REGION_SIZE=4000
    BASE_FILENAME=""
    BAM=""
    CHROMHMM=""
    CHROMS_NUM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    BIN_SIZE=50
    SHAPES_COMPREHENSIVE=""
    TRAINING=""
    while getopts n:d:c:i:r:s:a:u: option; do
        case "${option}" in
            n) CELL_LINE=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            c) CHROMHMM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            s) SHAPES_COMPREHENSIVE=$(realpath $OPTARG);;
            a) TRAINING=$(realpath $OPTARG);;
            u) CUTOFFS=$(realpath $OPTARG);;
        esac
    done
    BASE_PATH=$BASE_FILENAME/$CELL_LINE
    
#Print message to user.
    echo "Annotating with the following settings:"
    echo "Name is $CELL_LINE"
    echo "Base directory is $BASE_FILENAME"
    echo "Learning from $TRAINING"
    echo "Cutoffs in $CUTOFFS"
    echo "ChromHMM file is $CHROMHMM"
    echo "Bin size in WIG file is $BIN_SIZE"
    echo "Region size for annotation is $REGION_SIZE"
    echo "Shapes will be saved in $SHAPES_COMPREHENSIVE"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"
	
    for 
    
    #Create all needed directories.
    SOM_OUT="$BASE_PATH/som_output"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    SOM_OUT_FILTERED="$BASE_PATH/som_output_filtered"
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED="$BASE_PATH/som_output_shifted"
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL="$BASE_PATH/som_output_final"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG="$BASE_PATH/wig_chroms"
    SHAPE_ANNOTATED="$BASE_PATH/anno_beds"
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED="$BASE_PATH/anno_beds_sorted"
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL="$BASE_PATH/anno_beds_final"
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    if [[ ! -e ${SHAPES_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}_${CELL_LINE}
    fi
    SHAPES="$BASE_PATH/shapes"
    SHAPES_LOG="$BASE_PATH/shapes_log"
	
# Method for running the pipeline for a chromosome.
	run_pipeline() {
		local c=$1
		
		#Run the SOM.
		python vnssom.py $TRAINING/chrom$c.pkl $SOM_OUT/$c.pkl $CUTOFFS $TRAINING $REGION_SIZE $BIN_SIZE
		echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
		
		#Merge shifted regions.
		python merge_shifted.py $SOM_OUT/$c.pkl $SOM_OUT_SHIFTED/$c.pkl
		echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
		
		#Annotate regions with shape.
		python make_shape_bed.py $TRAINING/chrom$c.pkl $SOM_OUT_SHIFTED/$c.pkl $SHAPE_ANNOTATED/$c.bed
        echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"

        # python consolidate.py $SHAPE_ANNOTATED_SORTED/anno${c} $SHAPE_ANNOTATED_FINAL/anno${c}
        # cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/anno${c}.bed
        # awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/clusters_anno${c}
        # cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/scores_anno${c}.bed
        # echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
        
        #Save shapes to file.
        bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/$c.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/$c.bed			
        bedtools sort -i $CHROMHMM_INTERSECTS/$c.bed > $CHROMHMM_INTERSECTS/${c}_sorted.bed
        python find_chromhmm_distrib.py $CHROMHMM_INTERSECTS/${c}_sorted.bed $SOM_OUT_SHIFTED/$c.pkl $CHROMHMM_DISTRIB/$c.pkl
	}
    
    #Run the pipeline from split WIG files to final set of shapes.
    pids=""
    for f in $CHROMS_NUM;
        do 
            run_pipeline $f & 
            pids="$pids $!"
        done
        
    #Merge shapes entries across chromosomes.
    wait $pids
    python ../common_scripts/merge_significant.py $SHAPES_COMPREHENSIVE $SHAPES $SHAPES_LOG
        
    #Exit
	wait
	echo Done!
	exit 0