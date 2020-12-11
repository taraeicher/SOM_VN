#!/bin/bash
#SBATCH --job-name=learnShapesBCellHigh
#SBATCH --nodes=1 --ntasks=1 --mem=3g
#SBATCH --partition=gpu
#SBATCH --gres=gpu:k20x:1,lscratch:10
#SBATCH --time=1:00:00
#SBATCH --output=shapesBCellHigh_%A_%a.out
#SBATCH --error=shapesBCellHigh_%A_%a.err
#SBATCH --array=1-1000
#Variables
	BASE_PATH="/data/eichertd/som_vn_data"
	LSCRATCH_PATH="/lscratch/$SLURM_JOB_ID"
    CELL_LINE="b_cell_high"
    REGION_SIZE=4000
    CHROMHMM=$BASE_PATH/chromhmm/b_cell_E032.bed
    BIN_SIZE=50
	CHROM_IDX=$((SLURM_ARRAY_TASK_ID/100+1))
	ITER=$((SLURM_ARRAY_TASK_ID%100+1))

# Read chromosome using index.
	CHROM=$(head -$CHROM_IDX $BASE_PATH/chroms_rand_$ITER | tail -1)
	echo $CHROM

# Load modules.
	module load cuDNN/7.6.5/CUDA-10.1 python/3.7
	module load bedtools
    
#Create all needed directories.
    SOM_OUT=$LSCRATCH_PATH/som_${CELL_LINE}_${ITER}
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
    TRAINING_FILES=$BASE_PATH/training_${CELL_LINE}_shifted
	TRAINING_ANNOTATION_FILES=$BASE_PATH/training_annotation_${CELL_LINE}
    SOM_OUT_FILTERED=$BASE_PATH/som_${CELL_LINE}_${ITER}_filtered
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
    SOM_OUT_SHIFTED=$LSCRATCH_PATH/som_${CELL_LINE}_${ITER}_shifted
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
    SOM_OUT_FINAL=$LSCRATCH_PATH/som_${CELL_LINE}_${ITER}_final
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
    WIG=$BASE_PATH/wig/${CELL_LINE}
    SHAPE_ANNOTATED=$BASE_PATH/anno_${CELL_LINE}_${ITER}
    if [[ ! -e $SHAPE_ANNOTATED ]]; then
        mkdir $SHAPE_ANNOTATED
    fi
    SHAPE_ANNOTATED_SORTED=$LSCRATCH_PATH/anno_${CELL_LINE}_${ITER}_sorted
    if [[ ! -e $SHAPE_ANNOTATED_SORTED ]]; then
        mkdir $SHAPE_ANNOTATED_SORTED
    fi
    SHAPE_ANNOTATED_FINAL=$LSCRATCH_PATH/anno_${CELL_LINE}_${ITER}_final
    if [[ ! -e $SHAPE_ANNOTATED_FINAL ]]; then
        mkdir $SHAPE_ANNOTATED_FINAL
    fi
    CHROMHMM_INTERSECTS=$BASE_PATH/anno_${CELL_LINE}_${ITER}_intersect
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        mkdir $CHROMHMM_INTERSECTS
    fi
	SHAPES_COMPREHENSIVE=$BASE_PATH/anno_${CELL_LINE}_${ITER}_consolidated
    if [[ ! -e ${SHAPES_COMPREHENSIVE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}
    fi

	#Run the SOM.
	python som_vn.py $TRAINING_FILES/chrom${CHROM} $SOM_OUT/chrom${CHROM} $WIG/$CELL_LINE.chr${CHROM}.wig $REGION_SIZE $BIN_SIZE 0 False
	echo -e "---------------------------------------------SOM model is ready for chrom $c.-----------------------------------------\n"
	
	#Remove all shapes to which no regions map.
	python remove_by_cutoff.py $SOM_OUT/chrom${CHROM}som_centroid 1 $SOM_OUT_FILTERED/chrom${CHROM}som_centroid
	echo -e "------------------------------------------------Removal complete for chrom $c.---------------------------------------\n"
	
	#Merge shifted regions.
	python merge_shifted.py $SOM_OUT_FILTERED/chrom${CHROM}som_centroid $SOM_OUT_SHIFTED/chrom${CHROM}som_centroid 0
	echo -e "------------------------------------------------Merging complete for chrom $c.----------------------------------------\n"
	
	#Remove duplicate shapes using kmeans.
	python kmeans_shapes.py $SOM_OUT_SHIFTED/chrom${CHROM}som_centroid $SOM_OUT_FINAL/chrom${CHROM}som_centroid
	echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------\n"
	
	#Annotate regions with shape.
	python make_shape_bed.py $TRAINING_ANNOTATION_FILES/chrom${CHROM} $SOM_OUT_FINAL/chrom${CHROM}som_centroid $SHAPE_ANNOTATED/anno${CHROM} 0
	echo -e "\n------------------------------------Initial annotations complete for chrom $c.-----------------------------\n"
	
	bedtools sort -i  $SHAPE_ANNOTATED/anno${CHROM} > $SHAPE_ANNOTATED_SORTED/anno${CHROM}
	python consolidate.py $SHAPE_ANNOTATED_SORTED/anno${CHROM} $SHAPE_ANNOTATED_FINAL/anno${CHROM}
	cut -d$'\t' -f 1,2,3,4,5 $SHAPE_ANNOTATED_FINAL/anno${CHROM} > $SHAPE_ANNOTATED_FINAL/anno${CHROM}.bed
	awk '{ print $6}' $SHAPE_ANNOTATED_FINAL/anno${CHROM} > $SHAPE_ANNOTATED_FINAL/clusters_anno${CHROM}
	cut -d$'\t' -f 7,8,9,10 $SHAPE_ANNOTATED_FINAL/anno${CHROM} > $SHAPE_ANNOTATED_FINAL/scores_anno${CHROM}.bed
	echo -e "\n------------------------------------Consolidating complete for chrom $c.-----------------------------\n"
	
	#Save shapes to file.
	bedtools intersect -wao -a $SHAPE_ANNOTATED_FINAL/anno${CHROM}.bed -b $CHROMHMM > $CHROMHMM_INTERSECTS/anno${CHROM}.bed			
	bedtools sort -i $CHROMHMM_INTERSECTS/anno${CHROM}.bed > $CHROMHMM_INTERSECTS/anno${CHROM}_sorted.bed
	python consolidate_chromHMM.py $CHROMHMM_INTERSECTS/anno${CHROM}_sorted.bed $SOM_OUT_FINAL/chrom${CHROM}som_centroid ${SHAPES_COMPREHENSIVE} ${WIG}/${CELL_LINE}.chr${CHROM}.wig ${CHROM} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${CHROM} 0
	