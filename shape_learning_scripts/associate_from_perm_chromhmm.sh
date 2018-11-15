#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   

    SRC=""
    DEST=""
    BASE_FILENAME=""
    WIG=""
    ANNOTATION_WIG=""
    CHROMHMM=""
    CHROMHMM_PERM=""
    while getopts s:d:b:w:a:c:p: option; do
        case "${option}" in
            s) SRC=$OPTARG;;
            d) DEST=$OPTARG;;
            b) BASE_FILENAME=$(realpath $OPTARG);;
            w) WIG=$(realpath $OPTARG);;
            a) ANNOTATION_WIG=$(realpath $OPTARG);;
            c) CHROMHMM=$(realpath $OPTARG);;
            p) CHROMHMM_PERM=$(realpath $OPTARG);;
        esac
    done
    
    BASE_SRC="$BASE_FILENAME/$SRC"
    BASE_DEST="$BASE_FILENAME/$DEST"
    ANNOTATED_SHAPE="$BASE_SRC/anno_beds_final"
    TRAINING_ANNOTATION_FILES="$BASE_SRC/training_anno"
    ISECTS="$BASE_SRC/anno_intersects_chromhmm_perm"
    SOM="$BASE_SRC/som_output_final"
    SHAPE_LIST_ALL="$BASE_SRC/chromhmm_perm_percentage"
    SHAPE_LIST="$BASE_SRC/chromhmm_perm_merged"
    LOG="$BASE_SRC/ig_chromHmm_log_chromhmm_perm"
    ANNOTATED_CONSOLIDATED="$BASE_DEST/annotated_consolidated_perm_$SRC"
    ANNOTATED="$BASE_DEST/annotated_perm_$SRC"
    ANNOTATED_SORTED="$BASE_DEST/annotated_sorted_perm_$SRC"
    TO_ANNOTATE="$BASE_DEST/annotation_files"
    CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    
#Print message to user.
    echo "Annotating with the following settings:"
    echo "Source is $SRC"
    echo "Target is $DEST"
    echo "Base directory is $BASE_FILENAME"
    echo "ChromHMM file is $CHROMHMM"
    echo "Wig file is $WIG"
    echo "Wig file to annotate is $ANNOTATION_WIG"
    echo "Permuted ChromHMM file will be saved in in $CHROMHMM_PERM"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"

#Create all needed directories.
if [[ ! -e $ISECTS ]]; then
	mkdir $ISECTS
fi
if [[ ! -e $ANNOTATED ]]; then
	mkdir $ANNOTATED
fi
if [[ ! -e $ANNOTATED_SORTED ]]; then
	mkdir $ANNOTATED_SORTED
fi
if [[ ! -e $ANNOTATED_CONSOLIDATED ]]; then
	mkdir $ANNOTATED_CONSOLIDATED
fi
if [[ ! -e ${SHAPE_LIST_ALL}_${SRC} ]]; then
	mkdir ${SHAPE_LIST_ALL}_${SRC}
fi

#Build a shape list on a set of chromosomes.
#Loop through each chromosome.
for chr in $CHROMS;
	do 
        #Permute ChromHMM data.
        python permute_chromhmm.py $CHROMHMM $CHROMHMM_PERM
        
        #Find fake associations.
        bedtools intersect -wao -a $ANNOTATED_SHAPE/anno${chr}.bed -b $CHROMHMM_PERM > $ISECTS/anno${chr}_all.bed			
        bedtools sort -i $ISECTS/anno${chr}_all.bed > $ISECTS/anno${chr}_sorted_all.bed
        python consolidate_chromHMM.py $ISECTS/anno${chr}_sorted_all.bed $SOM/chrom${chr}som_centroid $SHAPE_LIST_ALL  ${WIG}${chr}.wig ${chr} $SRC $TRAINING_ANNOTATION_FILES/chrom${chr}window3 0
	done

python ../common_scripts/merge_significant.py $SHAPE_LIST_ALL $SHAPE_LIST $LOG
    
#Annotate the cell line and evaluate annotations.
annotate() {
        local c=$1
        #Annotate, sort, and consolidate.
        python ../annotation_scripts/make_annotated_bed.py $TO_ANNOTATE/chrom${c} $SHAPE_LIST $ANNOTATED/anno${c} $ANNOTATION_WIG/${DEST}.chr${c}.wig 0.0
        bedtools sort -i  $ANNOTATED/anno${c} > $ANNOTATED_SORTED/anno${c}
        bedtools sort -i  $ANNOTATED/anno${c}clust > $ANNOTATED_SORTED/anno${c}.clust
        python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED/anno${c} $ANNOTATED_CONSOLIDATED/anno${c}
        echo "Consolidated chrom $c"

        #Make BED file.
        cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/anno${c}.bed
        awk '{ print $6}' $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/clusters_anno${c}
        cut -d$'\t' -f 7,8,9,10 $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/scores_anno${c}.bed
        echo "Made bed for chrom $c"
    }
for f in $CHROMS;
    do 
        annotate $f & 
    done 
    
wait
echo Done!
exit 0