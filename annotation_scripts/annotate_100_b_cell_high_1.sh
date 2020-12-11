#!/bin/bash
#SBATCH --job-name=annotateBCellHigh
#SBATCH --nodes=1 --ntasks=1 --mem=3g
#SBATCH --gres=glscratch:10
#SBATCH --time=1:00:00
#SBATCH --output=annotateBCellHigh_%A_%a.out
#SBATCH --error=annotateBCellHigh_%A_%a.err
#SBATCH --array=1-1000
#Variables
	BASE_PATH="/data/eichertd/som_vn_data"
	LSCRATCH_PATH="/lscratch/$SLURM_JOB_ID"
    CELL_LINE="b_cell_high"
    REGION_SIZE=4000
    CHROMHMM=$BASE_PATH/chromhmm/b_cell_E032.bed
    BIN_SIZE=50
	CHROM_IDX=$((SLURM_ARRAY_TASK_ID/100))
	ITER=$((SLURM_ARRAY_TASK_ID%100))
	
# Read chromosome using index.
	CHROM=$(head -$CHROM_IDX $BASE_PATH/notchroms_rand_$ITER | tail -1)
	echo $CHROM

# Load modules.
	module load python/3.7
	module load bedtools

cd /home/eichertd/som_vn_code/SOM_VN/annotation_scripts

BASE_PATH="/data/eichertd/som_vn_data"
CELL_LINE="b_cell_high"
TO_ANNOTATE=$BASE_PATH/training_annotation_${CELL_LINE}
SHAPES=$BASE_PATH/shapes_${CELL_LINE}_${ITER}
CHROMHMM=$BASE_PATH/chromhmm/b_cell_E032.bed
WIG=$BASE_PATH/wig/${CELL_LINE}
ANNOTATED_TGT=$LSCRATCH_PATH/annotated_${CELL_LINE}_${ITER}
    if [[ ! -e $ANNOTATED_TGT ]]; then
        mkdir $ANNOTATED_TGT
    fi
ANNOTATED_SORTED_TGT=$LSCRATCH_PATH/annotated_sorted_${CELL_LINE}_${ITER}
    if [[ ! -e $ANNOTATED_SORTED_TGT ]]; then
        mkdir $ANNOTATED_SORTED_TGT
    fi
ANNOTATED_CONSOLIDATED_TGT=$BASE_PATH/annotated_consolidated_${CELL_LINE}_${ITER}
    if [[ ! -e $ANNOTATED_CONSOLIDATED_TGT ]]; then
        mkdir $ANNOTATED_CONSOLIDATED_TGT
    fi
ANNO_MERGED=$BASE_PATH/annotated_merged_${CELL_LINE}_${ITER}
        if [[ ! -e $ANNO_MERGED ]]; then
            mkdir $ANNO_MERGED
        fi

#Annotating the regions.
python make_annotated_bed.py $TO_ANNOTATE/chrom${CHROM} $SHAPES $ANNOTATED_TGT/anno${CHROM} $WIG/$CELL_LINE.chr${CHROM}.wig  0.0
    
#Sorting the annotated regions.
bedtools sort -i  $ANNOTATED_TGT/anno${CHROM} > $ANNOTATED_SORTED_TGT/anno${CHROM}
bedtools sort -i  $ANNOTATED_TGT/anno${CHROM}clust > $ANNOTATED_SORTED_TGT/anno${CHROM}.clust

#Consolidating the sorted annotated regions.
python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED_TGT/anno${CHROM} $ANNOTATED_CONSOLIDATED_TGT/anno${CHROM}
    
# Creating BED, cluster, and score files.
cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${CHROM} > $ANNOTATED_CONSOLIDATED_TGT/anno${CHROM}.bed
awk '{ print $6}' $ANNOTATED_CONSOLIDATED_TGT/anno${CHROM} > $ANNOTATED_CONSOLIDATED_TGT/clusters_anno${CHROM}
    
#Intersect our annotations with the ground truth.
bedtools intersect -wao -a $ANNOTATED_CONSOLIDATED_TGT/anno${CHROM}.bed -b $CHROMHMM > $ANNO_MERGED/anno${CHROM}.bed
