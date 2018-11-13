#PBS -l nodes=1:ppn=4
#PBS -l walltime=4:00:00 
#!/bin/bash   

#Enable job control.
	set -m

#Variables
    CELL_LINE="A549"
	BASE_PATH="/fs/project/PAS0272/Tara/DNase_SOM/$CELL_LINE"
	CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
    WINDOW_INDEX=3
    WINDOW_SIZE=4000
	BIN_SIZE=50
    CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${CELL_LINE}_chromhmm_15_liftOver.bed"
	
#Create all needed directories.
    cd "/fs/project/PAS0272/Tara/DNase_SOM/scripts"
	TRAINING_SHIFTED="$BASE_PATH/training_shifted_perm_below50"
    if [[ ! -e $TRAINING_SHIFTED ]]; then
        mkdir $TRAINING_SHIFTED
    fi
    TRAINING_FILES="$BASE_PATH/training_perm_below50"
    if [[ ! -e $TRAINING_FILES ]]; then
        mkdir $TRAINING_FILES
    fi
    TRAINING_ANNOTATION_FILES="$BASE_PATH/training_anno"
    if [[ ! -e $TRAINING_ANNOTATION_FILES ]]; then
        mkdir $TRAINING_ANNOTATION_FILES
    fi
    SOM_OUT="$BASE_PATH/som_output_perm_below50"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    SOM_OUT_FILTERED="$BASE_PATH/som_output_filtered_perm_below50"
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED="$BASE_PATH/som_output_shifted_perm_below50"
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL="$BASE_PATH/som_output_final_perm_below50"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG_ORIG="$BASE_PATH/wig_chroms"
    WIG="$BASE_PATH/wig_chroms_perm"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    SHAPE_ANNOTATED="$BASE_PATH/anno_beds_perm_below50"
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED="$BASE_PATH/anno_beds_sorted_perm_below50"
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL="$BASE_PATH/anno_beds_final_perm_below50"
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS="$BASE_PATH/anno_intersects_perm_below50"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
    DATABASE_COMPREHENSIVE="$BASE_PATH/database_all_perm_below50"
    if [[ ! -e ${DATABASE_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${DATABASE_COMPREHENSIVE}_${CELL_LINE}
    fi
    DISTRIB_FIGS="$BASE_PATH/annotation_distribution_figs_perm_below50"
    if [[ ! -e $DISTRIB_FIGS ]]; then
        mkdir $DISTRIB_FIGS
    fi
    DATABASE="$BASE_PATH/database_perm_below50"
    DATABASE_REAL="$BASE_PATH/database_under50"
    DATABASE_LOG="$BASE_PATH/database_log_perm_below50"
    PVALS="$BASE_PATH/pvalues_perm_below50"
    if [[ ! -e $PVALS ]]; then
        mkdir $PVALS
    fi
    CROSSCORR_REAL="$BASE_PATH/crosscorr_below50"
    if [[ ! -e $CROSSCORR_REAL ]]; then
        mkdir $CROSSCORR_REAL
    fi
    CROSSCORR="$BASE_PATH/crosscorr_perm_below50"
    if [[ ! -e $CROSSCORR ]]; then
        mkdir $CROSSCORR
    fi
    CLUSTER_FIGS="$BASE_PATH/cluster_fig_perm_below50"
    if [[ ! -e $CLUSTER_FIGS ]]; then
        mkdir $CLUSTER_FIGS
    fi
    TO_ANNOTATE="$BASE_PATH/annotation_files"
    
    gcc -pthread -lm -o runGetData getFileData_rpkmBelow50.c
    module load python/2.7
	
#Method for running the pipeline for a chromosome.
	run_pipeline() {
		local c=$1
        
        #Permute WIG.
        #python permute_wig.py ${WIG_ORIG}/${CELL_LINE}.chr${c}.wig ${WIG}/${CELL_LINE}.chr${c}.wig
        
		#Generate input files for training and annotation.
		# ./runGetData $WIG/$CELL_LINE.chr$c.wig $BIN_SIZE 0 Y $c $TRAINING_FILES/chrom$c
        # ./runGetData $WIG/$CELL_LINE.chr$c.wig $BIN_SIZE 0 N $c $TRAINING_ANNOTATION_FILES/chrom$c
		# echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------\n"
		
		# Shuffle the input files.
		# cd $TRAINING_FILES
		# shuf chrom${c}window${WINDOW_INDEX} > shuf_chrom${c}
        # cd "/fs/project/PAS0272/Tara/DNase_SOM/scripts"
		# echo -e "--------------------------------------------Shuffling complete for chrom $c.------------------------------------------\n"
			
		# Shift the input to its best representation.
		# python shift_input.py $TRAINING_FILES/shuf_chrom$c $TRAINING_SHIFTED/chrom$c $BIN_SIZE $WINDOW_SIZE $WIG/$CELL_LINE.chr$c.wig
		# echo -e "----------------------------------------------Shifting complete for chrom $c.-----------------------------------------\n"
		
		#Run the SOM.
		# python som_auto.py $TRAINING_SHIFTED/chrom$c $SOM_OUT/chrom$c $WIG/$CELL_LINE.chr${c}.wig $WINDOW_SIZE $BIN_SIZE
		# echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
		
		# #Remove all SOM centroids to which no regions map.
		# python remove_by_cutoff.py $SOM_OUT/chrom${c}som_centroid 1 $SOM_OUT_FILTERED/chrom${c}som_centroid
		# echo -e "------------------------------------------------Removal complete for chrom $c.---------------------------------------\n"
		
		# #Merge shifted regions.
		# python merge_shifted.py $SOM_OUT_FILTERED/chrom${c}som_centroid $SOM_OUT_SHIFTED/chrom${c}som_centroid
		# echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
		
		# #Remove duplicate centroids using kmeans.
		# python kmeans_centroids.py $SOM_OUT_SHIFTED/chrom${c}som_centroid $SOM_OUT_FINAL/chrom${c}som_centroid
		# echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------\n"
        
        # #Print out the clusters.
        # python print_clusters.py $SOM_OUT_FINAL/chrom${c}som_centroid $CLUSTER_FIGS/chrom${c} ${c}
		
		#Annotate regions with cluster data.
		# python make_cluster_bed.py $TRAINING_ANNOTATION_FILES/chrom${c}window${WINDOW_INDEX} $SOM_OUT_FINAL/chrom${c}som_centroid $SHAPE_ANNOTATED/anno$c
        # echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
	
		# bedtools sort -i  $SHAPE_ANNOTATED/anno${c} > $SHAPE_ANNOTATED_SORTED/anno${c}
        # python consolidate_each_window.py $SHAPE_ANNOTATED_SORTED/anno${c} $SHAPE_ANNOTATED_FINAL/anno${c}
        # cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/anno${c}.bed
        # awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/clusters_anno${c}
        # cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${c} > $SHAPE_ANNOTATED_FINAL/scores_anno${c}.bed
        # echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
        
        # #Build comprehensive database.
        # bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${c}.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/anno${c}.bed			
        # bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
        # python consolidate_chromHMM_peakOnly.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $SOM_OUT/chrom${c}som_centroid $DATABASE_COMPREHENSIVE ${WIG}/${CELL_LINE}.chr${c}.wig $DISTRIB_FIGS ${c} $WINDOW_SIZE $PVALS $TRAINING_ANNOTATION_FILES/chrom${c}window$WINDOW_INDEX $CELL_LINE
        # python merge_significant.py $DATABASE_COMPREHENSIVE $DATABASE $DATABASE_LOG
        
        #Evaluate the distribution of cross-correlations for this and "real" clusters.
        # python make_annotated_bed_crosscorr.py $TO_ANNOTATE/chrom${c}window$WINDOW_INDEX $DATABASE_REAL $CROSSCORR_REAL/anno${c} $WIG_ORIG/${CELL_LINE}.chr${c}.wig
        # python make_annotated_bed_crosscorr.py $TO_ANNOTATE/chrom${c}window$WINDOW_INDEX $DATABASE $CROSSCORR/anno${c} $WIG_ORIG/${CELL_LINE}.chr${c}.wig
	}
	#Run the pipeline for each chromosome separately.
    for f in $CHROMS;
        do 
            run_pipeline $f & 
        done
    #Plot the cross-correlation distribution.
    module load python/3.5
    python plot_crosscorr_distrib.py $CROSSCORR_REAL/anno $CROSSCORR/anno $BASE_PATH/${CELL_LINE}_crosscorr_distrib
    
	wait
	echo Done!
	exit 0