#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

CELL="A549"
out="h1"
out_cap="H1"
CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${CELL}_chromhmm_15_liftOver.bed"
CHROMHMM_SPLIT="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${CELL}"
SCRIPTS="/fs/project/PAS0272/Tara/DNase_SOM/scripts"
SOM="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/som_output_final"
DATABASE_ALL="/fs/project/PAS0272/Tara/DNase_SOM/${out_cap}/database_all"
DATABASE="/fs/project/PAS0272/Tara/DNase_SOM/${out_cap}/database"
VAL="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/cluster_validity" 
BASE="/fs/project/PAS0272/Tara/DNase_SOM/$CELL"
WIG="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/wig_chroms/$CELL.chr"
WIG_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/Brain/wig_chroms/"
WIG_A549="/fs/project/PAS0272/Tara/DNase_SOM/A549/wig_chroms/"
WIG_H1="/fs/project/PAS0272/Tara/DNase_SOM/H1/wig_chroms/"
TSS_COVERAGE_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/Brain/tss_coverage_${out}"
TSS_COVERAGE_A549="/fs/project/PAS0272/Tara/DNase_SOM/A549/tss_coverage_${out}"
TSS_COVERAGE_H1="/fs/project/PAS0272/Tara/DNase_SOM/H1/tss_coverage_${out}"
ANNOTATED_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/Brain/annotated_consolidated_${out}"
ANNOTATED_A549="/fs/project/PAS0272/Tara/DNase_SOM/A549/annotated_consolidated_${out}"
ANNOTATED_H1="/fs/project/PAS0272/Tara/DNase_SOM/H1/annotated_consolidated_${out}"
SPLIT_BED="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/$CELL"
LOG="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/ig_chromHmm_log"
FIGS="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/annotation_figs/"
ANNOTATED="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/annotated_consolidated_${out}"
ANNOTATED_WITH_TSS="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/annotated_combined_${out}"
ANNOTATED_ALL_WIN="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/anno_beds_final"
ANNO_MERGED="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/annotated_merged_${out}"
ANNO_MERGED_NONE="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/annotated_merged_none"
ROC="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/roc_${out}"
PRECISION_RECALL="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/precision_recall_${out}"
UNKNOWNS="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/unknown_${out}"
GENOME="$HOME/hg38.chrom.sizes_sorted"
TSS_PROMOTER_SRC="/fs/project/PAS0272/Tara/DNase_SOM/Brain/tss_promoter"
TSS_PROMOTER="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/tss_promoter"
TSS_PROMOTER_CLUST="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/tss_promoter_clusters"
TSS_MERGED="/fs/project/PAS0272/Tara/DNase_SOM/Brain/tss_merged_allchrom"
TSS="/fs/project/PAS0272/Tara/DNase_SOM/TSS_hg38_first3cols.bed"
CLUSTER_PLOTS="/fs/project/PAS0272/Tara/DNase_SOM/${out_cap}/cluster_plots"
MISCLASSIFIED_PLOTS="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/misclassified_${out}"
WIG_DISTRIB="/fs/project/PAS0272/Tara/DNase_SOM/wig_distribs"
ANNO_MERGED_PERM="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/annotated_merged_perm_${out}"
ANNOTATED_PERM="$BASE/annotated_consolidated_perm_${out}"
WINDOW_INDEX=3
WIG_CHROMS="$BASE/wig_chroms_annotation"
ANNOTATION_DISTRIBS="/fs/project/PAS0272/Tara/DNase_SOM/${out_cap}/annotation_distribs"

#Create all needed directories.
ISECTS="/fs/project/PAS0272/Tara/DNase_SOM/$CELL/anno_intersects"
if [[ ! -e $ISECTS ]]; then
	mkdir $ISECTS
fi
if [[ ! -e $FIGS ]]; then
	mkdir $FIGS
fi
if [[ ! -e $ANNO_MERGED ]]; then
	mkdir $ANNO_MERGED
fi
if [[ ! -e $ANNO_MERGED_PERM ]]; then
	mkdir $ANNO_MERGED_PERM
fi
if [[ ! -e $ROC ]]; then
	mkdir $ROC
fi
if [[ ! -e $TSS_PROMOTER ]]; then
	mkdir $TSS_PROMOTER
fi
if [[ ! -e $TSS_PROMOTER_CLUST ]]; then
	mkdir $TSS_PROMOTER_CLUST
fi
if [[ ! -e $TSS_MERGED ]]; then
	mkdir $TSS_MERGED
fi
if [[ ! -e $TSS_COVERAGE_BRAIN ]]; then
	mkdir $TSS_COVERAGE_BRAIN
fi
if [[ ! -e $TSS_COVERAGE_A549 ]]; then
	mkdir $TSS_COVERAGE_A549
fi
if [[ ! -e $TSS_COVERAGE_H1 ]]; then
	mkdir $TSS_COVERAGE_H1
fi
if [[ ! -e $CLUSTER_PLOTS ]]; then
	mkdir $CLUSTER_PLOTS
fi

if [[ ! -e $UNKNOWNS ]]; then
	mkdir $UNKNOWNS
fi

if [[ ! -e $PRECISION_RECALL ]]; then
	mkdir $PRECISION_RECALL
fi

if [[ ! -e $MISCLASSIFIED_PLOTS ]]; then
	mkdir $MISCLASSIFIED_PLOTS
fi
if [[ ! -e $ANNOTATED_WITH_TSS ]]; then
	mkdir $ANNOTATED_WITH_TSS
fi
if [[ ! -e $ANNOTATION_DISTRIBS ]]; then
	mkdir $ANNOTATION_DISTRIBS
fi
if [[ ! -e $SPLIT_BED ]]; then
	mkdir $SPLIT_BED
fi

# #Predict promoters based on TSS.
# bedtools slop -l 1500 -r 1499 -i $TSS -g $GENOME | sort -V -k1,1 -k2,2 > $TSS_PROMOTER/promoters.bed
# bedtools complement -i $TSS_PROMOTER/promoters.bed -g $GENOME > $TSS_PROMOTER/not_promoters.bed
# bedtools makewindows -b $TSS_PROMOTER/not_promoters.bed -w 3000 > $TSS_PROMOTER/not_promoters_win.bed
# awk '{printf("%s\t%s\t%s\tPromoter\n", $1, $2, $3)}' $TSS_PROMOTER/promoters.bed > $TSS_PROMOTER/promoters_labeled.bed
# awk '{printf("%s\t%s\t%s\tNot_Promoter\n", $1, $2, $3)}' $TSS_PROMOTER/not_promoters_win.bed > $TSS_PROMOTER/not_promoters_win_labeled.bed
# cat $TSS_PROMOTER/promoters_labeled.bed $TSS_PROMOTER/not_promoters_win_labeled.bed | sortBed -i > $TSS_PROMOTER/tss_promoter_pred.bed
# cd $TSS_PROMOTER
# awk '{close(f);f=$1}{print > f".bed"}' $TSS_PROMOTER/tss_promoter_pred.bed
# cd $SCRIPTS
# find $TSS_PROMOTER/ -type f -name '*chrM*' -delete
# find $TSS_PROMOTER/ -type f -name '*chrEBV*' -delete
# find $TSS_PROMOTER/ -type f -name '*random*' -delete
# find $TSS_PROMOTER/ -type f -name '*chrUn*' -delete
# find $TSS_PROMOTER/ -type f -name '*GL*' -delete
# find $TSS_PROMOTER/ -type f -name '*KI*' -delete
# find $TSS_PROMOTER/ -type f -name '*JH*' -delete
# find $TSS_PROMOTER/ -type f -name '*KB*' -delete
# gcc -pthread -lm -o runGetRegions getBedRegions.c

# #Get the precision and recall for each chromosome, for:
# #1. Our predictions
# #2. TSS-based predictions
# #3. Our predictions combined with TSS
 #for chr in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y';
	 #do       
        # #Intersect TSS-based annotations with the ground truth.
        #./runGetRegions $WIG_CHROMS/$CELL.chr${chr}.wig $TSS_PROMOTER_SRC/chr${chr}.bed ${chr} $TSS_PROMOTER_CLUST/clusters${chr} $TSS_PROMOTER/final${chr}.bed
        #bedtools intersect -wao -a $TSS_PROMOTER/final${chr}.bed -b $CHROMHMM > $TSS_MERGED/anno${chr}.bed
        
        # #Intersect our annotations with the ground truth.
         #bedtools intersect -wao -a $ANNOTATED/anno${chr}.bed -b $CHROMHMM > $ANNO_MERGED/anno${chr}.bed

        #Combine shape-based and TSS-based annotations and intersect with the ground truth.
        #python combine_prediction_beds.py $ANNOTATED/anno${chr}.bed $TSS_PROMOTER/final${chr}.bed $ANNOTATED_WITH_TSS/${chr}.bed
        #gcc -pthread -lm -o runGetLines getBedLines.c
        #./runGetLines $WIG_CHROMS/$CELL.chr${chr}.wig ${ANNOTATED_WITH_TSS}/${chr}.bed ${chr} ${ANNOTATED_WITH_TSS}/clusters${chr} ${ANNOTATED_WITH_TSS}/final${chr}.bed
        #bedtools intersect -wao -a $ANNOTATED_WITH_TSS/${chr}.bed -b $CHROMHMM > ${ANNO_MERGED}/anno${chr}_tss.bed
        
        # #Intersect the permuted data with the ground truth.
        #bedtools intersect -wao -a $ANNOTATED_PERM/anno${chr}.bed -b $CHROMHMM > ${ANNO_MERGED_PERM}/anno${chr}.bed
     #done
#Plot the precision and recall over all groups: ours, tss-based, combined, and permuted.