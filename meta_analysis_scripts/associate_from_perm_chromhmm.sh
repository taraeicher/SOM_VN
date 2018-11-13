#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

SRC="Brain"
DEST="Brain"
ANNOTATED_SHAPE="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/anno_beds_final"
BASE_SRC="/fs/project/PAS0272/Tara/DNase_SOM/$SRC"
BASE_DEST="/fs/project/PAS0272/Tara/DNase_SOM/$DEST"
FIGS="$BASE_DEST/figs_chromhmm_perm/"
PVALS="$BASE_DEST/pvals_chromhmm_perm"
TRAINING_ANNOTATION_FILES="$BASE_SRC/training_anno"
SCRIPTS="/fs/project/PAS0272/Tara/DNase_SOM/scripts"
ISECTS="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/anno_intersects_chromhmm_perm"
SOM="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/som_output_final"
DATABASE_ALL="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/chromhmm_perm_percentage"
DATABASE="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/chromhmm_perm_merged"
WIG="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/wig_chroms/$SRC.chr"
ANNOTATION_WIG="$BASE_DEST/wig_chroms_annotation"
LOG="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/ig_chromHmm_log_chromhmm_perm"
ANNOTATED_CONSOLIDATED="$BASE_DEST/annotated_consolidated_perm_$SRC"
ANNOTATED="$BASE_DEST/annotated_perm_brain"
ANNOTATED_SORTED="$BASE_DEST/annotated_sorted_perm_$SRC"
TO_ANNOTATE="$BASE_DEST/annotation_files"
WINDOW_INDEX=3
CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
CHROMS="X Y"

#Create all needed directories.
if [[ ! -e $ISECTS ]]; then
	mkdir $ISECTS
fi
if [[ ! -e $PVALS ]]; then
	mkdir $PVALS
fi
if [[ ! -e $FIGS ]]; then
	mkdir $FIGS
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
if [[ ! -e ${DATABASE_ALL}_${SRC} ]]; then
	mkdir ${DATABASE_ALL}_${SRC}
fi

#Build a database on a subset of chromosomes.
#Loop through each chromosome and each window size.
# for chr in $CHROMS;
	# do 
        # #Permute ChromHMM data.
        # CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${SRC}/chr${chr}.bed"
        # CHROMHMM_PERM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${SRC}/chr${chr}_perm.bed"
        # python permute_chromhmm.py $CHROMHMM $CHROMHMM_PERM
        
        # #Find fake associations.
        # bedtools intersect -wao -a $ANNOTATED_SHAPE/anno${chr}.bed -b $CHROMHMM_PERM > $ISECTS/anno${chr}_all.bed			
        # bedtools sort -i $ISECTS/anno${chr}_all.bed > $ISECTS/anno${chr}_sorted_all.bed
        # python consolidate_chromHMM_peakOnly.py $ISECTS/anno${chr}_sorted_all.bed $SOM/chrom${chr}som_centroid $DATABASE_ALL  ${WIG}${chr}.wig $FIGS ${chr} $WINDOW_INDEX $PVALS/${chr} $TRAINING_ANNOTATION_FILES/chrom${chr}window$WINDOW_INDEX $SRC
	# done

# python merge_significant.py $DATABASE_ALL $DATABASE $LOG
    
#Annotate the cell line and evaluate annotations.
annotate() {
        local c=$1
        #Annotate, sort, and consolidate.
        python make_annotated_bed.py $TO_ANNOTATE/chrom${c}window$WINDOW_INDEX $DATABASE $ANNOTATED/anno${c} $ANNOTATION_WIG/${DEST}.chr${c}.wig 0.0
        bedtools sort -i  $ANNOTATED/anno${c} > $ANNOTATED_SORTED/anno${c}
        bedtools sort -i  $ANNOTATED/anno${c}clust > $ANNOTATED_SORTED/anno${c}.clust
        python consolidate_bed.py $ANNOTATED_SORTED/anno${c} $ANNOTATED_CONSOLIDATED/anno${c}
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