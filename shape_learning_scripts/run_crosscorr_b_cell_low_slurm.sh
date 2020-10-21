#!/bin/bash
#SBATCH --job-name=runCrosscorrBCellLow
#SBATCH --nodes=1 --ntasks=1 --mem=1g --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=crosscorrBCellLow_%A_%a.out
#SBATCH --error=crosscorrBCellLow_%A_%a.err
#SBATCH --array=1-22

# Load modules.
	module load python/3.7
	
#Variables
	BASE_PATH="/data/eichertd/som_vn_data"
    CELL_LINE="b_cell_low"
	
#Create all needed directories.
    TRAINING_ANNOTATION_FILES=$BASE_PATH/training_annotation_${CELL_LINE}
    WIG=$BASE_PATH/wig/${CELL_LINE}
	WIG_ORIG=$WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}.wig
	WIG_PERM=$WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}_perm.wig
	SHAPES_REAL=$BASE_PATH/shapes_${CELL_LINE}
	SHAPES_PERM=$BASE_PATH/shapes_${CELL_LINE}_perm
	CROSSCORR_REAL=$BASE_PATH/crosscorr_${CELL_LINE}
    if [[ ! -e ${CROSSCORR_REAL} ]]; then
        mkdir ${CROSSCORR_REAL}
    fi
	CROSSCORR_PERM=$BASE_PATH/crosscorr_${CELL_LINE}_perm
    if [[ ! -e ${CROSSCORR_PERM} ]]; then
        mkdir ${CROSSCORR_PERM}
    fi
	
	# Run crosscorr.
	python ../annotation_scripts/make_annotated_bed_crosscorr.py $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} $SHAPES_PERM $CROSSCORR_PERM/anno${SLURM_ARRAY_TASK_ID} $WIG_ORIG 0
    python ../annotation_scripts/make_annotated_bed_crosscorr.py $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} $SHAPES_REAL $CROSSCORR_REAL/anno${SLURM_ARRAY_TASK_ID} $WIG_PERM 0
    echo -e "\n------------------------------------Finished annotation with crosscorr for chrom $c.-----------------------------\n"
	