#PBS -l nodes=1:ppn=28
#PBS -l walltime=4:00:00 
#!/bin/bash   

#Enable job control.
	set -m

#Variables
    CELL_LINE=""
	BASE_FILENAME=""
	CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
    REGION_SIZE=4000
	BIN_SIZE=50
    CHROMHMM=""
    WIG_ORIG=""
    while getopts n:d:c:i:r:w: option; do
        case "${option}" in
            n) CELL_LINE=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            c) CHROMHMM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            r) REGION_SIZE=$OPTARG;;
            w) WIG_ORIG=$(realpath $OPTARG);;
        esac
    done
    BASE_PATH=$BASE_FILENAME/$CELL_LINE
    PYTHON_VERSION=2.7	
	
#Create all needed directories.
	TRAINING_SHIFTED="$BASE_PATH/training_shifted_perm"
    if [[ ! -e $TRAINING_SHIFTED ]]; then
        mkdir $TRAINING_SHIFTED
    fi
    TRAINING_FILES="$BASE_PATH/training_perm"
    if [[ ! -e $TRAINING_FILES ]]; then
        mkdir $TRAINING_FILES
    fi
    TRAINING_ANNOTATION_FILES="$BASE_PATH/training_perm"
    if [[ ! -e $TRAINING_ANNOTATION_FILES ]]; then
        mkdir $TRAINING_ANNOTATION_FILES
    fi
    SOM_OUT="$BASE_PATH/som_output_perm"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    SOM_OUT_FILTERED="$BASE_PATH/som_output_filtered_perm"
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED="$BASE_PATH/som_output_shifted_perm"
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL="$BASE_PATH/som_output_final_perm"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG="$BASE_PATH/wig_chroms_perm"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    SHAPE_ANNOTATED="$BASE_PATH/anno_beds_perm"
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED="$BASE_PATH/anno_beds_sorted_perm"
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL="$BASE_PATH/anno_beds_final_perm"
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects_perm"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    SHAPES_COMPREHENSIVE="$BASE_PATH/SHAPES_all_perm"
    if [[ ! -e ${SHAPES_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}_${CELL_LINE}
    fi
    SHAPES="$BASE_PATH/shapes_perm"
    SHAPES_REAL="$BASE_PATH/shapes"
    SHAPES_LOG="$BASE_PATH/shapes_log_perm"
    PVALS="$BASE_PATH/pvalues_perm"
    CROSSCORR_REAL="$BASE_PATH/crosscorr"
    if [[ ! -e $CROSSCORR_REAL ]]; then
        mkdir $CROSSCORR_REAL
    fi
    CROSSCORR="$BASE_PATH/crosscorr_perm"
    if [[ ! -e $CROSSCORR ]]; then
        mkdir $CROSSCORR
    fi
    SHAPE_FIGS="$BASE_PATH/shape_fig_perm"
    if [[ ! -e $SHAPE_FIGS ]]; then
        mkdir $SHAPE_FIGS
    fi
    TO_ANNOTATE="$BASE_PATH/annotation_files"
    
    module load python/2.7
    gcc -pthread -lm -o run_get_data ../common_scripts/get_file_data.c
    echo -e "------------------------------------------------------Compilation complete.-------------------------------------------------\n"
	
#Method for running the pipeline for a chromosome.
	run_pipeline() {
		local c=$1
        
        #Permute WIG.
        python permute_wig.py ${WIG_ORIG}/${CELL_LINE}.chr${c}.wig ${WIG}/${CELL_LINE}.chr${c}.wig
        
		#Generate input files for training and annotation.
        ./run_get_data $WIG/$CELL_LINE.chr$c.wig $BIN_SIZE 0 Y $c $TRAINING_FILES/chrom${c} $REGION_SIZE
        ./run_get_data $WIG/$CELL_LINE.chr$c.wig $BIN_SIZE 0 Y $c $TRAINING_FILES/chrom${c} $REGION_SIZE
		echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------\n"
		
		#Shuffle the input files.
		shuf $TRAINING_FILES/chrom${c}window3 > $TRAINING_FILES/shuf_chrom${c}
		echo -e "--------------------------------------------Shuffling complete for chrom $c.------------------------------------------\n"
			
		#Shift the input to its best representation.
		python shift_input.py $TRAINING_FILES/shuf_chrom${c} $TRAINING_SHIFTED/chrom${c} $BIN_SIZE $REGION_SIZE $WIG/$CELL_LINE.chr$c.wig false 0
		echo -e "----------------------------------------------Shifting complete for chrom $c.-----------------------------------------\n"
		
		#Run the SOM.
        python som_vn.py $TRAINING_SHIFTED/chrom$c $SOM_OUT/chrom${c} $WIG/$CELL_LINE.chr$c.wig $REGION_SIZE $BIN_SIZE 0 False
		echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
		
		#Remove all SOM centroids to which no regions map.
		python remove_by_cutoff.py $SOM_OUT/chrom${c}som_centroid 1 $SOM_OUT_FILTERED/chrom${c}som_centroid
		echo -e "------------------------------------------------Removal complete for chrom $c.---------------------------------------\n"
		
		#Merge shifted regions.
		python merge_shifted.py $SOM_OUT_FILTERED/chrom${c}som_centroid $SOM_OUT_SHIFTED/chrom${c}som_centroid 0
		echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
		
		#Remove duplicate centroids using kmeans.
		python kmeans_shapes.py $SOM_OUT_SHIFTED/chrom${c}som_centroid $SOM_OUT_FINAL/chrom${c}som_centroid
		echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------\n"
		
		#Annotate regions with shape data.
		python make_shape_bed.py $TRAINING_ANNOTATION_FILES/chrom${c}window3 $SOM_OUT_FINAL/chrom${c}som_centroid $SHAPE_ANNOTATED/anno${c}test 0
        echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
	
		bedtools sort -i  $SHAPE_ANNOTATED/anno${c} > $SHAPE_ANNOTATED_SORTED/anno${c}
        python consolidate.py $SHAPE_ANNOTATED_SORTED/anno${c} $SHAPE_ANNOTATED_FINAL/anno${c}
        cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/anno${c}.bed
        awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/shapes_anno${c}
        cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/scores_anno${c}.bed
        echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
        
        #Build comprehensive shapes.
        bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${c}.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/anno${c}.bed			
        bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
        python consolidate_chromHMM.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $SOM_OUT_FINAL/chrom${c}som_centroid ${SHAPES_COMPREHENSIVE} ${WIG}/${CELL_LINE}.chr${c}.wig ${c} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${c}window3 0
        echo -e "\n------------------------------------Shapes saved for chrom $c.-----------------------------\n"
	}
	#Run the pipeline for each chromosome separately.
     for f in $CHROMS;
        do 
            run_pipeline $f &
        done
        
    #Merge shapes entries across chromosomes.
    python ../common_scripts/merge_significant.py $SHAPES_COMPREHENSIVE $SHAPES $SHAPES_LOG
    echo -e "\n------------------------------------Merged shapes across chromosomes.-----------------------------\n"
    
    #Evaluate the distribution of cross-correlations for this and "real" shapes.
    for c in $CHROMS;
        do
            python ../annotation_scripts/make_annotated_bed_crosscorr.py $TO_ANNOTATE/chrom${c} $SHAPES $CROSSCORR/anno${c} $WIG_ORIG/${CELL_LINE}.chr${c}.wig 0
            python ../annotation_scripts/make_annotated_bed_crosscorr.py $TO_ANNOTATE/chrom${c} $SHAPES_REAL $CROSSCORR_REAL/anno${c} $WIG_ORIG/${CELL_LINE}.chr${c}.wig 0
            echo -e "\n------------------------------------Finished annotation with crosscorr for chrom $c.-----------------------------\n"
        done
    
    #Plot the cross-correlation distribution.
    module load python/3.5
    python ../meta_analysis_scripts/plot_crosscorr_distrib.py $CROSSCORR_REAL/anno $CROSSCORR/anno $BASE_PATH/${CELL_LINE}_crosscorr_distrib ${CELL_LINE}
    echo -e "\n------------------------------------Plotted crosscorr distributions.-----------------------------\n"
    
	wait
	echo Done!
	exit 0