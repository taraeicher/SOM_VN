#PBS -l nodes=1:ppn=28
#PBS -l walltime=10:00:00
#!/bin/bash   
c=$1
ITER=$2
BASE=$3
CELL_LINE=$4

TO_ANNOTATE="$BASE/annotation_files" 
DATABASE="$BASE/database_${ITER}"
CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm"
WINDOW_INDEX=3
WIG="$BASE/wig_chroms"
ANNOTATED_TGT="$BASE/annotated_${ITER}"
    if [[ ! -e $ANNOTATED_TGT ]]; then
        mkdir $ANNOTATED_TGT
    fi
ANNOTATED_SORTED_TGT="$BASE/annotated_sorted_${ITER}"
    if [[ ! -e $ANNOTATED_SORTED_TGT ]]; then
        mkdir $ANNOTATED_SORTED_TGT
    fi
ANNOTATED_CONSOLIDATED_TGT="$BASE/annotated_consolidated_${ITER}"
    if [[ ! -e $ANNOTATED_CONSOLIDATED_TGT ]]; then
        mkdir $ANNOTATED_CONSOLIDATED_TGT
    fi
ANNO_MERGED="$BASE/annotated_merged_${ITER}"
        if [[ ! -e $ANNO_MERGED ]]; then
            mkdir $ANNO_MERGED
        fi
    
#Annotating the regions.
python make_annotated_bed.py $TO_ANNOTATE/chrom${c}window${WINDOW_INDEX} $DATABASE $ANNOTATED_TGT/anno${c} $WIG/${CELL_LINE}.chr${c}.wig 0.0
    
#Sorting the annotated regions.
bedtools sort -i  $ANNOTATED_TGT/anno${c} > $ANNOTATED_SORTED_TGT/anno${c}
bedtools sort -i  $ANNOTATED_TGT/anno${c}clust > $ANNOTATED_SORTED_TGT/anno${c}.clust

#Consolidating the sorted annotated regions.
python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED_TGT/anno${c} $ANNOTATED_CONSOLIDATED_TGT/anno${c}
    
# Creating BED, cluster, and score files.
cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${c} > $ANNOTATED_CONSOLIDATED_TGT/anno${c}.bed
cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED_TGT/anno${c}.clust > $ANNOTATED_CONSOLIDATED_TGT/anno${c}clust.bed
awk '{ print $6}' $ANNOTATED_CONSOLIDATED_TGT/anno${c} > $ANNOTATED_CONSOLIDATED_TGT/clusters_anno${c}
cut -d$'\t' -f 7,8,9,10 $ANNOTATED_CONSOLIDATED_TGT/anno${c} > $ANNOTATED_CONSOLIDATED_TGT/scores_anno${c}.bed
    
#Intersect our annotations with the ground truth.
bedtools intersect -wao -a $ANNOTATED_CONSOLIDATED_TGT/anno${c}.bed -b $CHROMHMM/${CELL_LINE}_chromhmm_15_liftOver.bed > $ANNO_MERGED/anno${c}.bed
