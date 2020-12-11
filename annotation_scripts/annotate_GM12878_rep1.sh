#!/bin/bash
#SBATCH --job-name=annotateGM12878Rep1
#SBATCH --nodes=1 --ntasks=1 --mem=3g
#SBATCH --time=2:00:00
#SBATCH --output=annotateGM12878_rep1_%A_%a.out
#SBATCH --error=annotateGM12878_rep1_%A_%a.err
#SBATCH --array=1-22

module load python/3.7
module load bedtools

BASE_PATH="/data/eichertd/som_vn_data"
CELL_LINE="GM12878_rep1"
TO_ANNOTATE=$BASE_PATH/training_annotation_${CELL_LINE}
SHAPES=$BASE_PATH/shapes_${CELL_LINE}
CHROMHMM=$BASE_PATH/chromhmm/GM12878_E116.bed
WIG=$BASE_PATH/wig/${CELL_LINE}
ANNOTATED_TGT="$BASE_PATH/annotated_${CELL_LINE}"
    if [[ ! -e $ANNOTATED_TGT ]]; then
        mkdir $ANNOTATED_TGT
    fi
ANNOTATED_SORTED_TGT="$BASE_PATH/annotated_sorted_${CELL_LINE}"
    if [[ ! -e $ANNOTATED_SORTED_TGT ]]; then
        mkdir $ANNOTATED_SORTED_TGT
    fi
ANNOTATED_CONSOLIDATED_TGT="$BASE_PATH/annotated_consolidated_${CELL_LINE}"
    if [[ ! -e $ANNOTATED_CONSOLIDATED_TGT ]]; then
        mkdir $ANNOTATED_CONSOLIDATED_TGT
    fi
ANNO_MERGED="$BASE_PATH/annotated_merged_${CELL_LINE}"
        if [[ ! -e $ANNO_MERGED ]]; then
            mkdir $ANNO_MERGED
        fi
    
#Annotating the regions.
python make_annotated_bed.py $TO_ANNOTATE/chrom${SLURM_ARRAY_TASK_ID} $SHAPES $ANNOTATED_TGT/anno${SLURM_ARRAY_TASK_ID} $WIG/$CELL_LINE.chr${SLURM_ARRAY_TASK_ID}.wig  0.0
    
#Sorting the annotated regions.
bedtools sort -i  $ANNOTATED_TGT/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_SORTED_TGT/anno${SLURM_ARRAY_TASK_ID}
bedtools sort -i  $ANNOTATED_TGT/anno${SLURM_ARRAY_TASK_ID}clust > $ANNOTATED_SORTED_TGT/anno${SLURM_ARRAY_TASK_ID}.clust

#Consolidating the sorted annotated regions.
python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED_TGT/anno${SLURM_ARRAY_TASK_ID} $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID}
    
# Creating BED, cluster, and score files.
cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID}.bed
cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID}.clust > $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID}clust.bed
awk '{ print $6}' $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_CONSOLIDATED_TGT/clusters_anno${SLURM_ARRAY_TASK_ID}

cut -d$'\t' -f 7,8,9,10 $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_CONSOLIDATED_TGT/scores_anno${SLURM_ARRAY_TASK_ID}.bed
    
#Intersect our annotations with the ground truth.
bedtools intersect -wao -a $ANNOTATED_CONSOLIDATED_TGT/anno${SLURM_ARRAY_TASK_ID}.bed -b $CHROMHMM > $ANNO_MERGED/anno${SLURM_ARRAY_TASK_ID}.bed
