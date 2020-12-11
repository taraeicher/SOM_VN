#!/bin/bash
#SBATCH --job-name=magnitudeBCellLow
#SBATCH --nodes=1 --ntasks=1 --mem=3g
#SBATCH --time=1:00:00
#SBATCH --output=magnitudeBCellLow_%A_%a.out
#SBATCH --error=magnitudeBCellLow_%A_%a.err
#SBATCH --array=1-22

# Load modules.
	module load cuDNN/7.6.5/CUDA-10.1 python/3.7
	module load bedtools

#Variables
	BASE_PATH="/data/eichertd/som_vn_data"
    CELL_LINE="b_cell_low"
    REGION_SIZE=4000
    CHROMHMM=$BASE_PATH/chromhmm/b_cell_E032.bed
    BIN_SIZE=50
    
#Create all needed directories.
    SOM_OUT=$BASE_PATH/som_${CELL_LINE}
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    TRAINING_FILES=$BASE_PATH/training_${CELL_LINE}_shifted
    WIG=$BASE_PATH/wig/${CELL_LINE}
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
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
	python make_magnitude_bed.py $TRAINING_FILES/chrom${SLURM_ARRAY_TASK_ID} $MAGNITUDE_ANNOTATED/anno${SLURM_ARRAY_TASK_ID} 0
	echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
	
	bedtools sort -i  $MAGNITUDE_ANNOTATED/anno${SLURM_ARRAY_TASK_ID} > $MAGNITUDE_ANNOTATED_SORTED/anno${SLURM_ARRAY_TASK_ID}.bed
	echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
	
	#Save shapes to file.
	bedtools intersect -wao -a $MAGNITUDE_ANNOTATED_SORTED/anno${SLURM_ARRAY_TASK_ID}.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}.bed			
	bedtools sort -i $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}.bed > $CHROMHMM_INTERSECTS/anno${SLURM_ARRAY_TASK_ID}_sorted.bed
	