#PBS -l nodes=1:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash  
#Enable job control.
	set -m

#Variables
    CELL_LINE=""
    BASE_FILENAME=""
    SHAPES=""
    BAM=""
    WIG_SPLIT_PATH=""
    BIN_SIZE=50
    REGION_SIZE=4000
    while getopts n:d:s: option; do
        case "${option}" in
            n) CELL_LINE=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            s) SHAPES=$(realpath $OPTARG);;
        esac
    done
	
    CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
	PYTHON_VERSION=2.7	
    WIG="$BASE_FILENAME/$CELL_LINE/$CELL_LINE.wig"
    TO_ANNOTATE="$BASE_FILENAME/$CELL_LINE/training_anno"

#Print message to user.
    echo "Annotating with the following settings:"
    echo "Name is $CELL_LINE"
    echo "Base directory is $BASE_FILENAME"
    echo "Shape file is $SHAPES"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"

#Create all needed directories.
    
	if [[ ! -e $TO_ANNOTATE ]]; then
		mkdir $TO_ANNOTATE
	fi
    WIGS="$BASE_FILENAME/$CELL_LINE/wig_chroms_annotation"
	if [[ ! -e $WIGS ]]; then
		mkdir WIGS
	fi
    ANNOTATED="$BASE_FILENAME/$CELL_LINE/annotated_nopromoter"
    if [[ ! -e $ANNOTATED ]]; then
		mkdir $ANNOTATED
	fi
    ANNOTATED_SORTED="$BASE_FILENAME/$CELL_LINE/annotated_sorted_nopromoter"
    if [[ ! -e $ANNOTATED_SORTED ]]; then
		mkdir $ANNOTATED_SORTED
	fi
    ANNOTATED_CONSOLIDATED="$BASE_FILENAME/$CELL_LINE/annotated_consolidated_nopromoter"
	if [[ ! -e $ANNOTATED_CONSOLIDATED ]]; then
		mkdir $ANNOTATED_CONSOLIDATED
	fi
    
    run_pipeline() {
        local c=$1
        
        #Annotating the regions.
        python make_annotated_bed.py $TO_ANNOTATE/chrom${c}_nopromoter_final $SHAPES $ANNOTATED/anno${c} $WIGS/${CELL_LINE}.chr${c}.wig 0.0
        echo -e "------------------------------------------------Annotation complete for chrom $c.------------------------------------------------\n"
        
        # Sorting the annotated regions.
        bedtools sort -i  $ANNOTATED/anno${c} > $ANNOTATED_SORTED/anno${c}
        bedtools sort -i  $ANNOTATED/anno${c}clust > $ANNOTATED_SORTED/anno${c}.clust

        # Consolidating the sorted annotated regions.
		python ../common_scripts/consolidate_bed.py $ANNOTATED_SORTED/anno${c} $ANNOTATED_CONSOLIDATED/anno${c}
		echo -e "------------------------------------------------Consolidation complete for chrom $c.----------------------------------------------\n"
		
        # Creating BED, cluster, and score files.
        cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/anno${c}.bed
        cut -d$'\t' -f 1,2,3,4,5 $ANNOTATED_CONSOLIDATED/anno${c}.clust > $ANNOTATED_CONSOLIDATED/anno${c}clust.bed
        awk '{ print $6}' $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/clusters_anno${c}
        cut -d$'\t' -f 7,8,9,10 $ANNOTATED_CONSOLIDATED/anno${c} > $ANNOTATED_CONSOLIDATED/scores_anno${c}.bed
    }
#Run the pipeline for each chromosome separately.
    for f in $CHROMS;
        do 
            run_pipeline $f & 
        done

	wait
	echo Done!
	exit 0
