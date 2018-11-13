#PBS -l nodes=1:ppn=28
#PBS -l walltime=10:00:00
#!/bin/bash   
#Directories used in analysis
CELL_LINE="H1"
BASE="/fs/project/PAS0272/Tara/DNase_SOM/${CELL_LINE}"
SCRIPTS="/fs/project/PAS0272/Tara/DNase_SOM/scripts"
cd $SCRIPTS

#Choose set of chromosomes randomly.
# chr=()
# for file in $(ls $BASE/som_output_final_${ITER}/); do
    # tmp=$(cut -d "s" -f 1 <<< "$file")
    # f=$(cut -d "m" -f 2 <<< "$tmp")
    # chr+=($f)
# done
# CHROMS=$(echo ${chr[*]})
CHROMS=$(echo "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" | tr " " "\n" | shuf -n 11)
NOTCHROMS=$(comm -3 <( echo $CHROMS | tr " " "\n" | sort ) <( echo "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" | tr " " "\n" | sort ))
WINDOW_INDEX=3
DATABASE_COMPREHENSIVE="$BASE/database_all_${ITER}"
    if [[ ! -e ${DATABASE_COMPREHENSIVE}_${CELL_LINE} ]]; then
        mkdir ${DATABASE_COMPREHENSIVE}_${CELL_LINE}
    fi
DATABASE="$BASE/database_${ITER}"
DATABASE_LOG="$BASE/database_log_${ITER}"
SOM_OUT="$BASE/som_output_${ITER}"
    if [[ ! -e $SOM_OUT ]]; then
        mkdir $SOM_OUT
    fi
SOM_OUT_FILTERED="$BASE/som_output_filtered_${ITER}"
    if [[ ! -e $SOM_OUT_FILTERED ]]; then
        mkdir $SOM_OUT_FILTERED
    fi
SOM_OUT_SHIFTED="$BASE/som_output_shifted_${ITER}"
    if [[ ! -e $SOM_OUT_SHIFTED ]]; then
        mkdir $SOM_OUT_SHIFTED
    fi
SOM_OUT_FINAL="$BASE/som_output_final_${ITER}"
    if [[ ! -e $SOM_OUT_FINAL ]]; then
        mkdir $SOM_OUT_FINAL
    fi
WIG="$BASE/wig_chroms"
    if [[ ! -e  $WIG ]]; then
        mkdir $WIG
    fi
ANNOTATED="$BASE/anno_beds_${ITER}"
    if [[ ! -e $ANNOTATED ]]; then
        mkdir $ANNOTATED
    fi
ANNOTATED_SORTED="$BASE/anno_beds_sorted_${ITER}"
    if [[ ! -e $ANNOTATED_SORTED ]]; then
        mkdir $ANNOTATED_SORTED
    fi
ANNOTATED_FINAL="$BASE/anno_beds_final_${ITER}"
    if [[ ! -e $ANNOTATED_FINAL ]]; then
        mkdir $ANNOTATED_FINAL
    fi
ANNOTATED_TGT="$BASE/annotated_consolidated_${ITER}"
    if [[ ! -e $ANNOTATED_TGT ]]; then
        mkdir $ANNOTATED_TGT
    fi
CHROMHMM_INTERSECTS="$BASE/anno_intersects_${ITER}"
    if [[ ! -e $CHROMHMM_INTERSECTS ]]; then
        echo "It's not here!"
        mkdir $CHROMHMM_INTERSECTS
    fi
ANNO_MERGED="$BASE/annotated_merged_${ITER}"
        if [[ ! -e $ANNO_MERGED ]]; then
            mkdir $ANNO_MERGED
        fi
PRECISION_RECALL="$BASE/precision_recall_${ITER}"
        if [[ ! -e $PRECISION_RECALL ]]; then
            mkdir $PRECISION_RECALL
        fi

#Run the pipeline for each chromosome separately.
cat /dev/null > $DATABASE_COMPREHENSIVE
for f in $CHROMS;
    do 
        ./learn_shapes.sh $f $ITER $BASE $CELL_LINE $SOM_OUT $SOM_OUT_FILTERED $SOM_OUT_SHIFTED $SOM_OUT_FINAL $WIG $ANNOTATED $ANNOTATED_SORTED $ANNOTATED_FINAL $CHROMHMM_INTERSECTS
        cat ${DATABASE_COMPREHENSIVE}_${f} >> $DATABASE_COMPREHENSIVE
    done
    
#Combine the data.
python merge_significant.py $DATABASE_COMPREHENSIVE $DATABASE $DATABASE_LOG

Annotate all remaning chromosomes.
for f in $NOTCHROMS;
    do 
        ./annotate.sh $f $ITER $BASE $CELL_LINE
    done
    
#Save the precision and recall for each chromosome.
python save_precision_recall.py $ANNO_MERGED/ $ANNOTATED_TGT/ $WIG/${CELL_LINE}.chr $PRECISION_RECALL/ $WINDOW_INDEX $CELL_LINE