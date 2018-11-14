#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   

CELL_LINE=""
GENOME=""
TSS=""
while getopts n:d:g:t: option; do
    case "${option}" in
        n) CELL_LINE=$OPTARG;;
        d) BASE_FILENAME=$(realpath $OPTARG);;
        g) GENOME=$(realpath $OPTARG);;
        t) TSS=$(realpath $OPTARG);;
    esac
done

#Create all needed directories.
TSS_PROMOTER="$BASE_FILENAME/$CELL_LINE/tss_promoter"
if [[ ! -e $TSS_PROMOTER ]]; then
	mkdir $TSS_PROMOTER
fi

#Predict promoters based on TSS.
bedtools slop -l 1500 -r 1499 -i $TSS -g $GENOME | sort -V -k1,1 -k2,2 > $TSS_PROMOTER/promoters.bed
bedtools complement -i $TSS_PROMOTER/promoters.bed -g $GENOME > $TSS_PROMOTER/not_promoters.bed
bedtools makewindows -b $TSS_PROMOTER/not_promoters.bed -w 3000 > $TSS_PROMOTER/not_promoters_win.bed
awk '{printf("%s\t%s\t%s\tPromoter\n", $1, $2, $3)}' $TSS_PROMOTER/promoters.bed > $TSS_PROMOTER/promoters_labeled.bed
awk '{printf("%s\t%s\t%s\tNot_Promoter\n", $1, $2, $3)}' $TSS_PROMOTER/not_promoters_win.bed > $TSS_PROMOTER/not_promoters_win_labeled.bed
cat $TSS_PROMOTER/promoters_labeled.bed $TSS_PROMOTER/not_promoters_win_labeled.bed | sortBed -i > $TSS_PROMOTER/tss_promoter_pred.bed
cd $TSS_PROMOTER
awk '{close(f);f=$1}{print > f".bed"}' $TSS_PROMOTER/tss_promoter_pred.bed
find $TSS_PROMOTER/ -type f -name '*chrM*' -delete
find $TSS_PROMOTER/ -type f -name '*chrEBV*' -delete
find $TSS_PROMOTER/ -type f -name '*random*' -delete
find $TSS_PROMOTER/ -type f -name '*chrUn*' -delete
find $TSS_PROMOTER/ -type f -name '*GL*' -delete
find $TSS_PROMOTER/ -type f -name '*KI*' -delete
find $TSS_PROMOTER/ -type f -name '*JH*' -delete
find $TSS_PROMOTER/ -type f -name '*KB*' -delete