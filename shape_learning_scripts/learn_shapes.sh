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
	
#Create all needed directories.
    SOM_OUT="$BASE_PATH/som_output_allchrom"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    TRAINING_FILES="$BASE_PATH/training_allchrom"
    if [[ ! -e $TRAINING_FILES ]]; then
        mkdir $TRAINING_FILES
    fi
    TRAINING_ANNOTATION_FILES="$BASE_PATH/training_anno_allchrom"
    if [[ ! -e $TRAINING_ANNOTATION_FILES ]]; then
        mkdir $TRAINING_ANNOTATION_FILES
    fi
    TRAINING_SHIFTED="$BASE_PATH/training_shifted_allchrom"
    if [[ ! -e $TRAINING_SHIFTED ]]; then
        mkdir $TRAINING_SHIFTED
    fi
    SOM_OUT_FILTERED="$BASE_PATH/som_output_filtered_allchrom"
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED="$BASE_PATH/som_output_shifted_allchrom"
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL="$BASE_PATH/som_output_final_allchrom"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG="$BASE_PATH/wig_chroms"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    SHAPE_ANNOTATED="$BASE_PATH/anno_beds_allchrom"
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED="$BASE_PATH/anno_beds_sorted_allchrom"
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL="$BASE_PATH/anno_beds_final_allchrom"
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects_allchrom"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    DISTRIB_FIGS="$BASE_PATH/annotation_distribution_figs_allchrom"
    if [[ ! -e $DISTRIB_FIGS ]]; then
        mkdir $DISTRIB_FIGS
    fi
    PVALS="$BASE_PATH/pvalues_allchrom"
    if [[ ! -e $PVALS ]]; then
        mkdir $PVALS
    fi
    if [[ ! -e ${SHAPES_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}_${CELL_LINE}
    fi
    SHAPES="$BASE_PATH/shapes_allchrom"
    SHAPES_LOG="$BASE_PATH/shapes_log_allchrom"

#Output RPKM intensities in a WIG file.
	gosr binbam -f 0 -n 1000 -t $CELL_LINE $BAM $BIN_SIZE $CELL_LINE > $BASE_PATH/$CELL_LINE.wig
	echo -e "------------------------------------------------------WIG file complete.-------------------------------------------------\n"
	
#Split WIG files by chromosome.
	python $WIG_SPLIT_PATH/wig_split.py $BASE_PATH/$CELL_LINE.wig $WIG/$CELL_LINE
	if [[ -e $WIG/*chrM* ]]; then
		rm $$WIG/*chrM*
	fi
	if [[ -e $WIG/*random* ]]; then
		rm $WIG/*random*
	fi
	if [[ -e $WIG/*chrEBV* ]]; then
		rm $WIG/*chrEBV*
	fi
    if [[ -e $WIG/*chrUn* ]]; then
		rm $WIG/*chrUn*
	fi
	if [[ -e $WIG/*_g* ]]; then
		rm $WIG/*_g*
	fi
    if [[ -e $WIG/*K* ]]; then
		rm $WIG/*K*
	fi
	echo -e "----------------------------------------------------WIG file split into chromosomes.------------------------------------\n"
	
#Data preprocessing
    gcc -pthread -lm -o run_get_data ../common_scripts/get_file_data.c
	echo -e "-----------------------------------------------------Compiled processing code.-----------------------------------------\n"
	
#Method for running the pipeline for a chromosome.
	run_pipeline() {
		local c=$1

		#Generate input files for training and annotation.
        ./run_get_data $WIGS/$CELL_LINE.chr$c.wig $BIN_SIZE 0 Y $c $TRAINING_FILES/chrom${c} $REGION_SIZE
        ./run_get_data $WIGS/$CELL_LINE.chr$c.wig $BIN_SIZE 0 N $c $TRAINING_ANNOTATION_FILES/chrom${c} $REGION_SIZE
		echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------\n"
		
		#Shuffle the input files.
		shuf $TRAINING_FILES/chrom${c} > $TRAINING_FILES/shuf_chrom${c}
		echo -e "--------------------------------------------Shuffling complete for chrom $c.------------------------------------------\n"
		
		#Shift the input to its best representation.
		python shift_input.py $TRAINING_FILES/shuf_chrom$c $TRAINING_SHIFTED/chrom$c $BIN_SIZE $REGION_SIZE $WIG/$CELL_LINE.chr$c.wig
		echo -e "----------------------------------------------Shifting complete for chrom $c.-----------------------------------------\n"
		
		#Run the SOM.
		python som_vn.py $TRAINING_SHIFTED/chrom$c $SOM_OUT/chrom$c $WIG/$CELL_LINE.chr$c.wig $REGION_SIZE $BIN_SIZE
		echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
		
		#Remove all shapes to which no regions map.
		python remove_by_cutoff.py $SOM_OUT/chrom${c}som_centroid 1 $SOM_OUT_FILTERED/chrom${c}som_centroid
		echo -e "------------------------------------------------Removal complete for chrom $c.---------------------------------------\n"
		
		#Merge shifted regions.
		python merge_shifted.py $SOM_OUT_FILTERED/chrom${c}som_centroid $SOM_OUT_SHIFTED/chrom${c}som_centroid
		echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
		
		#Remove duplicate shapes using kmeans.
		python kmeans_shapes.py $SOM_OUT_SHIFTED/chrom${c}som_centroid $SOM_OUT_FINAL/chrom${c}som_centroid
		echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------\n"
		
		#Annotate regions with shape.
		python make_shape_bed.py $TRAINING_ANNOTATION_FILES/chrom${c} $SOM_OUT_FINAL/chrom${c}som_centroid $SHAPE_ANNOTATED/anno$c
        echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
        
        bedtools sort -i  $SHAPE_ANNOTATED/anno${c} > $SHAPE_ANNOTATED_SORTED/anno${c}
        python consolidate.py $SHAPE_ANNOTATED_SORTED/anno${c} $SHAPE_ANNOTATED_FINAL/anno${c}
        cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/anno${c}.bed
        awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/clusters_anno${c}
        cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/scores_anno${c}.bed
        echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
        
        #Save shapes to file.
        bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${c}.bed -b $CHROMHMM/${CELL_LINE}_chromhmm_15_liftOver.bed > $CHROMHMM_INTERSECTS/anno${c}.bed			
        bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
        python consolidate_chromHMM.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $SOM_OUT_FINAL/chrom${c}som_centroid $SHAPES_COMPREHENSIVE ${WIG}/${CELL_LINE}.chr${c}.wig $DISTRIB_FIGS ${c} $REGION_SIZE $PVALS $TRAINING_ANNOTATION_FILES/chrom${c}window${WINDOW_INDEX} $CELL_LINE
	}
    
    #Run the pipeline from split WIG files to final set of annotations.
    for f in $CHROMS_NUM;
        do 
            run_pipeline $f & 
        done
        
    #Merge SHAPES entries across chromosomes.
    python merge_significant.py $SHAPES_COMPREHENSIVE $SHAPES $SHAPES_LOG
        
    #Exit
	wait
	echo Done!
	exit 0