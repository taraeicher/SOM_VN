#Enable job control.
	set -m

#Variables
	CELL_LINE="Brain"
	BASE_FILENAME="/fs/project/PAS0272/Tara/DNase_SOM"
    CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
	WIG_SPLIT_PATH="$HOME/taolib/Scripts/"
	BIN_SIZE=50
	PYTHON_VERSION=2.7	
    WINDOW_INDEX=3
    TSS_DIR="/fs/project/PAS0272/Tara/DNase_SOM/$CELL_LINE/tss_promoter"
	
#Create all needed directories.
    BAM="$BASE_FILENAME/bam/$CELL_LINE.bam"
    WIG="$BASE_FILENAME/$CELL_LINE/$CELL_LINE.wig"
    TO_ANNOTATE="$BASE_FILENAME/$CELL_LINE/annotation_files"
    DATABASE="$BASE_FILENAME/Brain/database"
    
	if [[ ! -e $TO_ANNOTATE ]]; then
		mkdir $TO_ANNOTATE
	fi
    WIGS="$BASE_FILENAME/$CELL_LINE/wig_chroms_annotation"
	if [[ ! -e $WIGS ]]; then
		mkdir WIGS
	fi
    ANNOTATED="$BASE_FILENAME/$CELL_LINE/annotated_H1_tss_and"
    if [[ ! -e $ANNOTATED ]]; then
		mkdir $ANNOTATED
	fi
    ANNOTATED_SORTED="$BASE_FILENAME/$CELL_LINE/annotated_H1_sorted_tss_and"
    if [[ ! -e $ANNOTATED_SORTED ]]; then
		mkdir $ANNOTATED_SORTED
	fi
    ANNOTATED_CONSOLIDATED="$BASE_FILENAME/$CELL_LINE/annotated_H1_consolidated_tss_and"
	if [[ ! -e $ANNOTATED_CONSOLIDATED ]]; then
		mkdir $ANNOTATED_CONSOLIDATED
	fi
    
#Output RPKM intensities in a WIG file.
	# gosr binbam -f 0 -n 1000 -t $CELL_LINE $BAM $BIN_SIZE $CELL_LINE > $WIG
	# echo -e "------------------------------------------------------WIG file complete.-------------------------------------------------\n"

#Split WIG files by chromosome.
	# python $WIG_SPLIT_PATH/wig_split.py $WIG $WIGS/$CELL_LINE
	# find $WIGS/ -type f -name '*chrM*' -delete
	# find $WIGS/ -type f -name '*chrEBV*' -delete
	# find $WIGS/ -type f -name '*random*' -delete
	# find $WIGS/ -type f -name '*chrUn*' -delete
	# find $WIGS/ -type f -name '*GL*' -delete
	# find $WIGS/ -type f -name '*KI*' -delete
	# echo -e "-----------------------------------------------------Splitting complete.------------------------------------------------\n"
    
#Split TSS by chromosome.
# file_list=$(ls $TSS_DIR/final*)
# for f in $file_list; 
    # do
        # chrbed=$(echo $f | awk -F"final" '{print $2}')
        # awk '{if($4 == "Promoter"){print $1,$2,$3}}' $TSS_DIR/final$chrbed > $TSS_DIR/promoter$chrbed
    # done
# cd $BASE_FILENAME/scripts	
#Data preprocessing
	# gcc -pthread -lm -o runGetData getFileData.c
	# echo -e "-----------------------------------------------------Compiling complete.------------------------------------------------\n"

    run_pipeline() {
        local c=$1
        
        # #Formatting the regions for annotation.
        # ./runGetData $WIGS/$CELL_LINE.chr$c.wig $BIN_SIZE 0 N $c $TO_ANNOTATE/chrom${c}
        # echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------------\n"
        
        # #Annotating the regions.
        # python make_annotated_bed_intersectTSS.py $TO_ANNOTATE/chrom${c}window${WINDOW_INDEX} $DATABASE $ANNOTATED/anno${c} $WIGS/${CELL_LINE}.chr${c}.wig 0.0 ${TSS_DIR}/promoter${c}.bed
        # echo -e "------------------------------------------------Annotation complete for chrom $c.------------------------------------------------\n"
        
        #Sorting the annotated regions.
        bedtools sort -i  $ANNOTATED/anno${c} > $ANNOTATED_SORTED/anno${c}
        bedtools sort -i  $ANNOTATED/anno${c}clust > $ANNOTATED_SORTED/anno${c}.clust

        #Consolidating the sorted annotated regions.
		 python consolidate_bed.py $ANNOTATED_SORTED/anno${c} $ANNOTATED_CONSOLIDATED/anno${c}
		 echo -e "------------------------------------------------Consolidation complete for chrom $c.------------------------------------------------\n"
		
        # Creating BED, cluster, and score files.
        cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/anno${c}.bed
        cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c}.clust > $ANNOTATED_CONSOLIDATED/anno${c}clust.bed
        awk '{ print $6}' $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/clusters_anno${c}
        cut -d$'\t' -f 7,8,9,10 $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/scores_anno${c}.bed
		echo -e "\n------------------------------------Consolidating per window size complete for chrom $c.-----------------------------------------"
    }
#Run the pipeline for each chromosome separately.
    for f in $CHROMS;
        do 
            run_pipeline $f & 
        done

	wait
	echo Done!
	exit 0
