#!/bin/bash
#SBATCH --job-name=permChromHMMH1High
#SBATCH --nodes=1 --ntasks=1 --mem=1g --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=permChromHMMH1High_%A_%a.out
#SBATCH --error=permChromHMMH1High_%A_%a.err
#SBATCH --array=1-22
    
# Load modules.
	module load python/3.7
	module load bedtools
	
#Variables
	BASE_PATH="/data/eichertd/som_vn_data"
	CELL_LINE="H1_high"
	CHROMHMM_PERM=$BASE_PATH/chromhmm/H1_E003_perm.bed
	TRAINING_ANNOTATION_FILES=$BASE_PATH/training_annotation_${CELL_LINE}
	SHAPE_ANNOTATED_FINAL=$BASE_PATH/anno_${CELL_LINE}_final
	SOM_OUT_FINAL=$BASE_PATH/som_${CELL_LINE}_final
	WIG=$BASE_PATH/wig/${CELL_LINE}/${CELL_LINE}.chr${SLURM_ARRAY_TASK_ID}.wig

# Run
	#Create all needed directories.
	CHROMHMM_INTERSECTS=$BASE_PATH/anno_${CELL_LINE}_intersect_chromhmm_perm
	if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
		mkdir $CHROMHMM_INTERSECTS
	fi
	SHAPES_COMPREHENSIVE=$BASE_PATH/anno_${CELL_LINE}_consolidated_chromhmm_perm
	if [[ ! -e ${SHAPES_COMPREHENSIVE} ]]; then
		mkdir ${SHAPES_COMPREHENSIVE}
	fi

	#Find fake associations.
	bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}.bed -b $CHROMHMM_PERM > $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}.bed			
	bedtools sort -i $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}.bed > $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}_sorted.bed
	python consolidate_chromHMM.py $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}_sorted.bed $SOM_OUT_FINAL/chrom${SLURM_ARRAY_TASK_ID}som_centroid ${SHAPES_COMPREHENSIVE} $WIG ${SLURM_ARRAY_TASK_ID} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} 0

		