#!/bin/bash
#SBATCH --job-name=learnShapesGM12878Rep2
#SBATCH --nodes=1 --ntasks=1 --mem=3g
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k20x:1
#SBATCH --time=1:00:00
#SBATCH --output=shapesGM12878_Rep2_%A_%a.out
#SBATCH --error=shapesGM12878_Rep2_%A_%a.err
#SBATCH --array=1-22

# Load modules.
	module load cuDNN/7.6.5/CUDA-10.1 python/3.7
	module load bedtools

#Variables
	BASE_PATH="/data/eichertd/som_vn_data"
    CELL_LINE="GM12878_rep2"
    REGION_SIZE=4000
    CHROMHMM=$BASE_PATH/chromhmm/GM12878_E116.bed
    BIN_SIZE=50
    
#Create all needed directories.
    SOM_OUT=$BASE_PATH/som_${CELL_LINE}
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    TRAINING_FILES=$BASE_PATH/training_${CELL_LINE}
    if [[ ! -e $TRAINING_FILES ]]; then
        mkdir $TRAINING_FILES
    fi
    TRAINING_ANNOTATION_FILES=$BASE_PATH/training_annotation_${CELL_LINE}
    if [[ ! -e $TRAINING_ANNOTATION_FILES ]]; then
        mkdir $TRAINING_ANNOTATION_FILES
    fi
    TRAINING_SHIFTED=$BASE_PATH/training_${CELL_LINE}_shifted
    if [[ ! -e $TRAINING_SHIFTED ]]; then
        mkdir $TRAINING_SHIFTED
    fi
    SOM_OUT_FILTERED=$BASE_PATH/som_${CELL_LINE}_filtered
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED=$BASE_PATH/som_${CELL_LINE}_shifted
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL=$BASE_PATH/som_${CELL_LINE}_final
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG=$BASE_PATH/wig/${CELL_LINE}
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
    SHAPE_ANNOTATED=$BASE_PATH/anno_${CELL_LINE}
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED=$BASE_PATH/anno_${CELL_LINE}_sorted
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL=$BASE_PATH/anno_${CELL_LINE}_final
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS=$BASE_PATH/anno_${CELL_LINE}_intersect
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
	SHAPES_COMPREHENSIVE=$BASE_PATH/anno_${CELL_LINE}_consolidated
    if [[ ! -e ${SHAPES_COMPREHENSIVE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}
    fi


	#Generate input files for training and annotation.
	./run_get_data $WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}.wig $BIN_SIZE 0 Y ${SLURM_ARRAY_TASK_ID} $TRAINING_FILES/chrom${SLURM_ARRAY_TASK_ID} $REGION_SIZE
	./run_get_data $WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}.wig $BIN_SIZE 0 N ${SLURM_ARRAY_TASK_ID} $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} $REGION_SIZE
	echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------\n"
	
	#Shuffle the input files.
	shuf $TRAINING_FILES/chrom${SLURM_ARRAY_TASK_ID} > $TRAINING_FILES/shuf_chrom${SLURM_ARRAY_TASK_ID}
	echo -e "--------------------------------------------Shuffling complete for chrom $c.------------------------------------------\n"
	
	#Shift the input to its best representation.
	python shift_input.py $TRAINING_FILES/shuf_chrom${SLURM_ARRAY_TASK_ID} $TRAINING_SHIFTED/chrom${SLURM_ARRAY_TASK_ID} $BIN_SIZE $REGION_SIZE $WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}.wig false 0
	echo -e "----------------------------------------------Shifting complete for chrom $c.-----------------------------------------\n"
	
	#Run the SOM.
	python som_vn.py $TRAINING_SHIFTED/chrom${SLURM_ARRAY_TASK_ID} $SOM_OUT/chrom${SLURM_ARRAY_TASK_ID} $WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}.wig $REGION_SIZE $BIN_SIZE 0 False
	echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
	
	#Remove all shapes to which no regions map.
	python remove_by_cutoff.py $SOM_OUT/chrom${SLURM_ARRAY_TASK_ID}som_centroid 1 $SOM_OUT_FILTERED/chrom${SLURM_ARRAY_TASK_ID}som_centroid
	echo -e "------------------------------------------------Removal complete for chrom $c.---------------------------------------\n"
	
	#Merge shifted regions.
	python merge_shifted.py $SOM_OUT_FILTERED/chrom${SLURM_ARRAY_TASK_ID}som_centroid $SOM_OUT_SHIFTED/chrom${SLURM_ARRAY_TASK_ID}som_centroid 0
	echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
	
	#Remove duplicate shapes using kmeans.
	python kmeans_shapes.py $SOM_OUT_SHIFTED/chrom${SLURM_ARRAY_TASK_ID}som_centroid $SOM_OUT_FINAL/chrom${SLURM_ARRAY_TASK_ID}som_centroid
	echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------\n"
	
	#Annotate regions with shape.
	python make_shape_bed.py $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} $SOM_OUT_FINAL/chrom${SLURM_ARRAY_TASK_ID}som_centroid $SHAPE_ANNOTATED/anno${SLURM_ARRAY_TASK_ID} 0
	echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
	
	bedtools sort -i  $SHAPE_ANNOTATED/anno${SLURM_ARRAY_TASK_ID} > $SHAPE_ANNOTATED_SORTED/anno${SLURM_ARRAY_TASK_ID}
	python consolidate.py $SHAPE_ANNOTATED_SORTED/anno${SLURM_ARRAY_TASK_ID} $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}
	cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID} > $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}.bed
	awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID} > $SHAPE_ANNOTATED_FINAL/clusters_anno${SLURM_ARRAY_TASK_ID}
	cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID} > $SHAPE_ANNOTATED_FINAL/scores_anno${SLURM_ARRAY_TASK_ID}.bed
	echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
	
	#Save shapes to file.
	bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}.bed			
	bedtools sort -i $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}.bed > $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}_sorted.bed
	python consolidate_chromHMM.py $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}_sorted.bed $SOM_OUT_FINAL/chrom${SLURM_ARRAY_TASK_ID}som_centroid ${SHAPES_COMPREHENSIVE} ${WIG}/${CELL_LINE}.chr${SLURM_ARRAY_TASK_ID}.wig ${SLURM_ARRAY_TASK_ID} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} 0
	