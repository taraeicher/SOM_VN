#PBS -l nodes=2:ppn=28
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
    CELL_LINE="Brain"
    WINDOW_SIZE=4000
    WINDOW_INDEX=3
    BASE_PATH="/fs/project/PAS0272/Tara/DNase_SOM/$CELL_LINE"
    BAM="$BASE_PATH/$CELL_LINE/$CELL_LINE.bam"
    CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm"
    SCRIPTS="/fs/project/PAS0272/Tara/DNase_SOM/scripts"
    CHROMS_NUM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    WIG_SPLIT_PATH="$HOME/taolib/Scripts/"
    BIN_SIZE=50
	
#Create all needed directories.
    SOM_OUT="$BASE_PATH/som_output_allchrom_log"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    TRAINING_FILES="$BASE_PATH/training_allchrom_log"
    if [[ ! -e $TRAINING_FILES ]]; then
        mkdir $TRAINING_FILES
    fi
    TRAINING_ANNOTATION_FILES="$BASE_PATH/training_anno_allchrom_log"
    if [[ ! -e $TRAINING_ANNOTATION_FILES ]]; then
        mkdir $TRAINING_ANNOTATION_FILES
    fi
    TRAINING_SHIFTED="$BASE_PATH/training_shifted_allchrom_log"
    if [[ ! -e $TRAINING_SHIFTED ]]; then
        mkdir $TRAINING_SHIFTED
    fi
    SOM_OUT_FILTERED="$BASE_PATH/som_output_filtered_allchrom_log"
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED="$BASE_PATH/som_output_shifted_allchrom_log"
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL="$BASE_PATH/som_output_final_allchrom_log"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG_PRE="$BASE_PATH/wig_chroms"
    if [[ ! -e  $WIG_PRE ]]; then
        mkdir $WIG_PRE
    fi
    WIG="$BASE_PATH/wig_chroms_log"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    SHAPE_ANNOTATED="$BASE_PATH/anno_beds_allchrom_log"
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED="$BASE_PATH/anno_beds_sorted_allchrom_log"
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL="$BASE_PATH/anno_beds_final_allchrom_log"
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects_allchrom_log"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    DISTRIB_FIGS="$BASE_PATH/annotation_distribution_figs_allchrom_log"
    if [[ ! -e $DISTRIB_FIGS ]]; then
        mkdir $DISTRIB_FIGS
    fi
    PVALS="$BASE_PATH/pvalues_allchrom_log"
    if [[ ! -e $PVALS ]]; then
        mkdir $PVALS
    fi
    DATABASE_COMPREHENSIVE="$BASE_PATH/database_all_allchrom_log"
    if [[ ! -e ${DATABASE_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${DATABASE_COMPREHENSIVE}_${CELL_LINE}
    fi
    DATABASE="$BASE_PATH/database_allchrom_log"
    DATABASE_LOG="$BASE_PATH/database_log_allchrom_log"
    LIM=0.000005

# #Output RPKM intensities in a WIG file.
	# gosr binbam -f 0 -n 1000 -t $CELL_LINE $BAM $BIN_SIZE $CELL_LINE > $BASE_PATH/$CELL_LINE.wig
	# echo -e "------------------------------------------------------WIG file complete.-------------------------------------------------\n"
	
# # #Split WIG files by chromosome.
	# python $WIG_SPLIT_PATH/wig_split.py $BASE_PATH/$CELL_LINE.wig $WIG_PRE/$CELL_LINE
	# if [[ -e $WIG_PRE/*chrM* ]]; then
		# rm $$WIG_PRE/*chrM*
	# fi
	# if [[ -e $WIG_PRE/*random* ]]; then
		# rm $WIG_PRE/*random*
	# fi
	# if [[ -e $WIG_PRE/*chrEBV* ]]; then
		# rm $WIG_PRE/*chrEBV*
	# fi
    # if [[ -e $WIG_PRE/*chrUn* ]]; then
		# rm $WIG_PRE/*chrUn*
	# fi
	# if [[ -e $WIG_PRE/*_g* ]]; then
		# rm $WIG_PRE/*_g*
	# fi
    # if [[ -e $WIG_PRE/*K* ]]; then
		# rm $WIG_PRE/*K*
	# fi
	# echo -e "----------------------------------------------------WIG file split into chromosomes.------------------------------------\n"
	
#Data preprocessing
	#gcc -pthread -lm -o runGetData getFileData.c
	# echo -e "-----------------------------------------------------Compiled processing code.-----------------------------------------\n"
	
#Method for running the pipeline for a chromosome.
	run_pipeline() {
		local c=$1

        #Log-scale the WIG file.
        #python log_transform.py $WIG_PRE/$CELL_LINE.chr$c.wig $WIG/$CELL_LINE.chr$c.wig $LIM
        
		#Generate input files for training and annotation.
        # echo $WIG/$CELL_LINE.chr$c.wig        
		#./runGetData $WIG/$CELL_LINE.chr$c.wig $BIN_SIZE $(echo "l($LIM)" | bc -l) Y $c $TRAINING_FILES/chrom$c
		#./runGetData $WIG/$CELL_LINE.chr$c.wig $BIN_SIZE $(echo "l($LIM)" | bc -l) N $c $TRAINING_ANNOTATION_FILES/chrom$c
		#echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------\n"
		
		# #Shuffle the input files.
		# cd $TRAINING_FILES
		# shuf chrom${c}window${WINDOW_INDEX} > shuf_chrom${c}
        # cd $SCRIPTS
		# echo -e "--------------------------------------------Shuffling complete for chrom $c.------------------------------------------\n"
		
		#Shift the input to its best representation.
		# python shift_input.py $TRAINING_FILES/shuf_chrom$c $TRAINING_SHIFTED/chrom$c $BIN_SIZE $WINDOW_SIZE $WIG/$CELL_LINE.chr$c.wig 0.995 True $(echo "l($LIM)" | bc -l)
		# echo -e "----------------------------------------------Shifting complete for chrom $c.-----------------------------------------\n"
		
		#Run the SOM.
        # module load python/2.7
		# python som_auto.py $TRAINING_SHIFTED/chrom$c $SOM_OUT/chrom$c $WIG/$CELL_LINE.chr$c.wig $WINDOW_SIZE $BIN_SIZE $(echo "l($LIM)" | bc -l) 0.995 True
		# echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
		
		#Remove all SOM centroids to which no regions map.
		# python remove_by_cutoff.py $SOM_OUT/chrom${c}som_centroid 1 $SOM_OUT_FILTERED/chrom${c}som_centroid
		# echo -e "------------------------------------------------Removal complete for chrom $c.---------------------------------------\n"
		
		#Merge shifted regions.
		# python merge_shifted.py $SOM_OUT_FILTERED/chrom${c}som_centroid $SOM_OUT_SHIFTED/chrom${c}som_centroid $(echo "l($LIM)" | bc -l)
		# echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
		
		# #Remove duplicate centroids using kmeans.
		# python kmeans_centroids.py $SOM_OUT_SHIFTED/chrom${c}som_centroid $SOM_OUT_FINAL/chrom${c}som_centroid
		# echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------\n"
		
		# #Annotate regions with cluster data.
		# python make_cluster_bed.py $TRAINING_ANNOTATION_FILES/chrom${c}window${WINDOW_INDEX} $SOM_OUT_FINAL/chrom${c}som_centroid $SHAPE_ANNOTATED/anno$c $(echo "l($LIM)" | bc -l)
        # echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
        
        # bedtools sort -i  $SHAPE_ANNOTATED/anno${c} > $SHAPE_ANNOTATED_SORTED/anno${c}
        # python consolidate_each_window.py $SHAPE_ANNOTATED_SORTED/anno${c} $SHAPE_ANNOTATED_FINAL/anno${c}
        # cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/anno${c}.bed
        # awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/clusters_anno${c}
        # cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/scores_anno${c}.bed
        # echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
        
        #Build comprehensive database.
        # bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${c}.bed -b $CHROMHMM/${CELL_LINE}_chromhmm_15_liftOver.bed > $CHROMHMM_INTERSECTS/anno${c}.bed			
        # bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
        python consolidate_chromHMM_peakOnly.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $SOM_OUT_FINAL/chrom${c}som_centroid $DATABASE_COMPREHENSIVE ${WIG}/${CELL_LINE}.chr${c}.wig ${c} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${c}window${WINDOW_INDEX} $(echo "l($LIM)" | bc -l)
	}
    
    #Run the pipeline from split WIG files to final set of annotations.
    for f in $CHROMS_NUM;
        do 
            run_pipeline $f & 
        done
        
    #Merge database entries across chromosomes.
    #python merge_significant.py $DATABASE_COMPREHENSIVE $DATABASE $DATABASE_LOG
        
    #Exit
	wait
	echo Done!
	exit 0