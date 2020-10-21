#!/bin/bash
#SBATCH --job-name=cagtHeLaHigh
#SBATCH --nodes=1 --ntasks=1 --mem=1g --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=CAGT_HeLa_high_%A_%a.out
#SBATCH --error=CAGT_HeLa_high_%A_%a.err
#SBATCH --array=1-22

# Load modules.
	module load python/3.7
	#module load matlab
	module load bedtools
#Variables
	CELL_LINE="HeLa_high"
    BIN_SIZE=10
	BASE_PATH=/data/eichertd/som_vn_data/
    CAGT_PATH=/home/eichertd/som_vn_code/SOM_VN/cagt/trunk/matlab/src/
	CAGT_OUT=$BASE_PATH/cagt_out_${CELL_LINE}
	CHROMHMM=$BASE_PATH/chromhmm/hela_E117.bed
	TRAINING_SHIFTED=$BASE_PATH/training_${CELL_LINE}_shifted
	WIG=$BASE_PATH/wig/$CELL_LINE
	TRAINING_ANNOTATION_FILES=$BASE_PATH/training_annotation_${CELL_LINE}
    MERGE_DIST=0.8
    ITERATIONS=1000
    K=40
	
    #Move to the directory containing the scripts.
    cd /home/eichertd/som_vn_code/SOM_VN/shape_learning_scripts
    
    #Create all needed directories.

    MATLAB_MATRIX=$BASE_PATH/matlab_matrix_${CELL_LINE}
    if [[ ! -e $MATLAB_MATRIX ]]; then
        mkdir $MATLAB_MATRIX
    fi
    CAGT_OUT_SHIFTED=$BASE_PATH/cagt_out_shifted_${CELL_LINE}
    if [[ ! -e $CAGT_OUT_SHIFTED ]]; then
        mkdir $CAGT_OUT_SHIFTED
    fi
    ANNOTATED=$BASE_PATH/cagt_anno_beds_${CELL_LINE}
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
	ANNOTATED_SORTED=$BASE_PATH/cagt_anno_beds_sorted_${CELL_LINE}
    if [[ ! -e $ANNOTATED_SORTED ]]; then
        mkdir $ANNOTATED_SORTED
    fi
	ANNOTATED_FINAL=$BASE_PATH/cagt_anno_beds_final_${CELL_LINE}
    if [[ ! -e $ANNOTATED_FINAL ]]; then
        mkdir $ANNOTATED_FINAL
    fi
    INTERSECTS=$BASE_PATH/cagt_intersects_${CELL_LINE}
    if [[ ! -e $INTERSECTS ]]; then
        mkdir $INTERSECTS
    fi
    SHAPES_COMPREHENSIVE=$BASE_PATH/anno_${CELL_LINE}_consolidated_cagt
    if [[ ! -e ${SHAPES_COMPREHENSIVE} ]]; then
        mkdir ${SHAPES_COMPREHENSIVE}
    fi

    #Extract the signal and run CAGT.
    #python extract_signal.py $TRAINING $TRAINING_CSV/$CHROM.csv
    #module load matlab
    #matlab -nodisplay -nodesktop -r "run_cagt('$TRAINING_SHIFTED/chrom${SLURM_ARRAY_TASK_ID},'$MATLAB_MATRIX/${SLURM_ARRAY_TASK_ID}.mat','${SLURM_ARRAY_TASK_ID}','$CAGT_OUT/${SLURM_ARRAY_TASK_ID}.csv','$CAGT_PATH', '$MERGE_DIST', '$ITERATIONS', '$K')"
	#matlab -nodisplay -nodesktop -r "run_cagt('/data/eichertd/som_vn_data/training_b_cell_low_shifted/chrom1','/data/eichertd/som_vn_data/matlab_matrix_b_cell_low','1','/data/eichertd/som_vn_data/cagt_out_b_cell_low/1.csv','/home/eichertd/som_vn_code/SOM_VN/cagt/trunk/matlab/src/', '0.8', '1000', '40')"
    
    #Merge shifted regions.
	python merge_shifted.py $CAGT_OUT/${SLURM_ARRAY_TASK_ID}.csv $CAGT_OUT_SHIFTED/chrom${SLURM_ARRAY_TASK_ID} 0
    echo -e "Merging complete for chrom $CHROM.\n"
    
    #Annotate regions with shape.
	python make_shape_bed.py $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} $CAGT_OUT_SHIFTED/chrom${SLURM_ARRAY_TASK_ID}   $ANNOTATED/anno${SLURM_ARRAY_TASK_ID} 0
    echo -e "Initial annotations complete for chrom ${SLURM_ARRAY_TASK_ID}.\n"
	
	bedtools sort -i  $ANNOTATED/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_SORTED/anno${SLURM_ARRAY_TASK_ID}
	python consolidate.py $ANNOTATED_SORTED/anno${SLURM_ARRAY_TASK_ID} $ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}
	cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}.bed
	awk '{ print $6}' $ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_FINAL/clusters_anno${SLURM_ARRAY_TASK_ID}
	cut -d$'\t' -f 7,8,9,10 $ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID} > $ANNOTATED_FINAL/scores_anno${SLURM_ARRAY_TASK_ID}.bed
    
    #Intersect regions with ChromHMM.
    bedtools intersect -wao -a $ANNOTATED_FINAL/anno${SLURM_ARRAY_TASK_ID}.bed -b $CHROMHMM > $INTERSECTS/${SLURM_ARRAY_TASK_ID}.bed
    bedtools sort -i $INTERSECTS/${SLURM_ARRAY_TASK_ID}.bed > $INTERSECTS/${SLURM_ARRAY_TASK_ID}_sorted.bed
    python consolidate_chromHMM.py $INTERSECTS/${SLURM_ARRAY_TASK_ID}_sorted.bed $CAGT_OUT_SHIFTED/chrom${SLURM_ARRAY_TASK_ID} ${SHAPES_COMPREHENSIVE} ${WIG}/${CELL_LINE}.chr${SLURM_ARRAY_TASK_ID}.wig ${SLURM_ARRAY_TASK_ID} $CELL_LINE $TRAINING_ANNOTATION_FILES/chrom${SLURM_ARRAY_TASK_ID} 0
	
