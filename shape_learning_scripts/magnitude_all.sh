# Load modules.
	module load cuDNN/7.6.5/CUDA-10.1 python/3.7
	module load bedtools
	
CELL_LINES=(Brain_high H1_low A549 b_cell_low b_cell_high Brain_low H1_high HeLa_low HeLa_high heart_low heart_high stomach_low stomach_high GM12878)	
CHROMHMMS=(brain_E081.bed H1_E003.bed A549_E114.bed b_cell_E032.bed b_cell_E032.bed brain_E081.bed H1_E003.bed hela_E117.bed hela_E117.bed heart_E083.bed heart_E083.bed stomach_E092.bed stomach_E092.bed GM12878_E116.bed)
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
		TRAINING_FILES=$BASE_PATH/training_annotation_${CELL_LINE}
		WIG=$BASE_PATH/wig/${CELL_LINE}
		MAGNITUDE_ANNOTATED=$BASE_PATH/anno_${CELL_LINE}_magnitude
		if [[ ! -e $MAGNITUDE_ANNOTATED ]]; then
			mkdir $MAGNITUDE_ANNOTATED
		fi
		MAGNITUDE_ANNOTATED_SORTED=$BASE_PATH/anno_${CELL_LINE}_magnitude_sorted
		if [[ ! -e $MAGNITUDE_ANNOTATED_SORTED ]]; then
			mkdir $MAGNITUDE_ANNOTATED_SORTED
		fi
		CHROMHMM_INTERSECTS=$BASE_PATH/anno_${CELL_LINE}_magnitude_intersect
		if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
			mkdir $CHROMHMM_INTERSECTS
		fi
		MAGNITUDE_COMPREHENSIVE=$BASE_PATH/anno_${CELL_LINE}_magnitude_consolidated
		if [[ ! -e ${MAGNITUDE_COMPREHENSIVE} ]]; then
			mkdir ${MAGNITUDE_COMPREHENSIVE}
		fi

		
		#Annotate regions with magnitude.
		python make_magnitude_bed.py $TRAINING_FILES/chrom${c} $MAGNITUDE_ANNOTATED/anno${c} 0
		bedtools sort -i  $MAGNITUDE_ANNOTATED/anno${c} > $MAGNITUDE_ANNOTATED_SORTED/anno${c}.bed
		
		#Save shapes to file.
		bedtools intersect -wao -a $MAGNITUDE_ANNOTATED_SORTED/anno${c}.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/anno${c}.bed			
		bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
		
		#Consolidate.
		python consolidate_chromHMM_magnitude.py $CHROMHMM_INTERSECTS/anno${c}.bed  $MAGNITUDE_COMPREHENSIVE/anno${c} $WIG/${CELL_LINE}.chr${c}.wig ${c} ${CELL_LINE} ${TRAINING_FILES}/chrom${c} 0
	done 
done 
