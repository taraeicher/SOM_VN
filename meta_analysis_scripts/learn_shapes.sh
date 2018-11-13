#PBS -l nodes=1:ppn=28
#PBS -l walltime=20:00:00
#!/bin/bash   
c=$1
ITER=$2
BASE=$3
CELL_LINE=$4
SOM_OUT=$5
SOM_OUT_FILTERED=$6
SOM_OUT_SHIFTED=$7
SOM_OUT_FINAL=$8
WIG=$9
ANNOTATED=${10}
ANNOTATED_SORTED=${11}
ANNOTATED_FINAL=${12}
CHROMHMM_INTERSECTS=${13}

INPUT="$BASE/training"
INPUT_ANNO="$BASE/training_anno"
CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm"
WINDOW_SIZE=4000
WINDOW=3
BIN_SIZE=60
DATABASE="$BASE/database_all_${ITER}_${c}"
    if [[ ! -e ${DATABASE}_${CELL_LINE} ]]; then
        mkdir ${DATABASE}_${CELL_LINE}
    fi
module load python/2.7   
IMG="$BASE/consolidation_img_${ITER}"
    if [[ ! -e $IMG ]]; then
        mkdir $IMG
    fi
    
#Run the SOM.
python som_auto.py $INPUT/chrom${c}window${WINDOW} $SOM_OUT/chrom$c $WIG/$CELL_LINE.chr$c.wig $WINDOW_SIZE $BIN_SIZE

#Remove all SOM centroids to which no regions map.
python remove_by_cutoff.py $SOM_OUT/chrom${c}som_centroid 1 $SOM_OUT_FILTERED/chrom${c}som_centroid

#Merge shifted regions.
python merge_shifted.py $SOM_OUT_FILTERED/chrom${c}som_centroid $SOM_OUT_SHIFTED/chrom${c}som_centroid

#Remove duplicate centroids using kmeans.
python kmeans_centroids.py $SOM_OUT_SHIFTED/chrom${c}som_centroid $SOM_OUT_FINAL/chrom${c}som_centroid

#Annotate regions with cluster data.
python make_cluster_bed.py $INPUT_ANNO/chrom${c}window${WINDOW} $SOM_OUT_FINAL/chrom${c}som_centroid $ANNOTATED/anno$c

#Consolidate region annotations.
bedtools sort -i  $ANNOTATED/anno${c} > $ANNOTATED_SORTED/anno${c}
python consolidate_each_window.py $ANNOTATED_SORTED/anno${c} $ANNOTATED_FINAL/anno${c}
cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_FINAL/anno${c} > $ANNOTATED_FINAL/anno${c}.bed
awk '{ print $6}' $ANNOTATED_FINAL/anno${c} > $ANNOTATED_FINAL/clusters_anno${c}
cut -d$'\t' -f 7,8,9,10 $ANNOTATED_FINAL/anno${c} > $ANNOTATED_FINAL/scores_anno${c}.bed

# Build database for chromosome.
bedtools intersect -wao -a $ANNOTATED_FINAL/anno${c}.bed -b $CHROMHMM/${CELL_LINE}_chromhmm_15_liftOver.bed > $CHROMHMM_INTERSECTS/anno${c}.bed
bedtools sort -i $CHROMHMM_INTERSECTS/anno${c}.bed > $CHROMHMM_INTERSECTS/anno${c}_sorted.bed
python consolidate_chromHMM_peakOnly.py $CHROMHMM_INTERSECTS/anno${c}_sorted.bed $SOM_OUT_FINAL/chrom${c}som_centroid $DATABASE ${WIG}/$CELL_LINE.chr$c.wig ${c} $CELL_LINE ${INPUT_ANNO}/chrom${c}window${WINDOW}
