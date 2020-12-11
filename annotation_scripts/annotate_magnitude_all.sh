# Load modules.
	module load cuDNN/7.6.5/CUDA-10.1 python/3.7
	module load bedtools
	
CELL_LINES=(b_cell_low b_cell_high Brain_low H1_high HeLa_low HeLa_high heart_low heart_high stomach_low stomach_high GM12878)	
CHROMHMMS=(b_cell_E032.bed b_cell_E032.bed brain_E081.bed H1_E003.bed hela_E117.bed hela_E117.bed heart_E083.bed heart_E083.bed stomach_E092.bed stomach_E092.bed GM12878_E116.bed)
for i in {0..10};
do

#Variables
	 BASE_PATH="/data/eichertd/som_vn_data"
     CELL_LINE=${CELL_LINES[$i]}
     REGION_SIZE=4000
     CHROMHMM=$BASE_PATH/chromhmm/${CHROMHMMS[$i]}
     BIN_SIZE=50
	
	for c in {1..22}; do
    
	#Create all needed directories.
		 TO_ANNOTATE=$BASE_PATH/training_annotation_${CELL_LINE}
		 MAGNITUDES=$BASE_PATH/magnitudes_${CELL_LINE}
		 WIG=$BASE_PATH/wig/${CELL_LINE}
		 ANNOTATED_TGT=$BASE_PATH/annotated_magnitude_${CELL_LINE}
		# if [[ ! -e $ANNOTATED_TGT ]]; then
			# mkdir $ANNOTATED_TGT
		# fi
		 ANNOTATED_SORTED_TGT=$BASE_PATH/annotated_sorted_magnitude_${CELL_LINE}
		# if [[ ! -e $ANNOTATED_SORTED_TGT ]]; then
			# mkdir $ANNOTATED_SORTED_TGT
		# fi
		 ANNOTATED_CONSOLIDATED_TGT=$BASE_PATH/annotated_consolidated_magnitude_${CELL_LINE}
		# if [[ ! -e $ANNOTATED_CONSOLIDATED_TGT ]]; then
			# mkdir $ANNOTATED_CONSOLIDATED_TGT
		# fi
		 ANNO_MERGED=$BASE_PATH/annotated_merged_magnitude_${CELL_LINE}
		# if [[ ! -e ${ANNO_MERGED} ]]; then
			# mkdir ${ANNO_MERGED}
		# fi

		# #Annotate regions with magnitude.
		# python make_annotated_bed_magnitude.py $TO_ANNOTATE/chrom${c} $MAGNITUDES $ANNOTATED_TGT/anno${c} $WIG/$CELL_LINE.chr${c}.wig  0.0
		# bedtools sort -i  $ANNOTATED_TGT/anno${c} > $ANNOTATED_SORTED_TGT/anno${c}
		# bedtools sort -i  $ANNOTATED_TGT/anno${c}clust > $ANNOTATED_SORTED_TGT/anno${c}.clust
		# python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED_TGT/anno${c} $ANNOTATED_CONSOLIDATED_TGT/anno${c}
		
		# #Create bed and cluster files.
		# cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${c} > $ANNOTATED_CONSOLIDATED_TGT/anno${c}.bed
		# cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${c}.clust > $ANNOTATED_CONSOLIDATED_TGT/anno${c}clust.bed
		# awk '{print $6}' $ANNOTATED_CONSOLIDATED_TGT/anno${c} > $ANNOTATED_CONSOLIDATED_TGT/clusters_anno${c}
		
		#Intersect results with ground truth.
		bedtools intersect -wao -a $ANNOTATED_CONSOLIDATED_TGT/anno${c}.bed -b $CHROMHMM > $ANNO_MERGED/anno${c}.bed
	done 
done 