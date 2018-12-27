#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

#Directories that do not need to be created for each source-dest combination.
 SRC=(A549 A549 A549 Brain Brain Brain H1 H1 H1)
 DEST=(A549 Brain H1 A549 Brain H1 A549 Brain H1)
 TSS_PROMOTER_SRC="/fs/project/PAS0272/Tara/DNase_SOM/Brain/tss_promoter"
 
 #Compile the c code.
 gcc -pthread -lm -o run_get_regions ../common_scripts/get_bed_regions.c
 gcc -pthread -lm -o run_get_lines ../common_scripts/get_bed_lines.c
        
# Get the precision and recall for a subset of Brain tissue chromosomes for:
# 1. Our predictions
# 2. TSS-based predictions
# 3. Our predictions combined with TSS
# 4. Permuted predictions
# 5. Predictions using RPKM intensity

#Run analysis for each source and destination pair.
for i in 0 1 2 3 4 5 6 7 8;
    do
        src_cell=${SRC[$i]}
        dest_cell=${DEST[$i]}
        #Directories used in analysis
        BASE="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}"
        WIG_CHROMS="$BASE/wig_chroms_annotation"
        TSS_PROMOTER_SHAPES="$BASE/tss_promoter_shapes"
        if [[ ! -e $TSS_PROMOTER_SHAPES ]]; then
            mkdir $TSS_PROMOTER_SHAPES
        fi
        TSS_PROMOTER="$BASE/tss_promoter"
        if [[ ! -e $TSS_PROMOTER ]]; then
            mkdir $TSS_PROMOTER
        fi
        ANNOTATED_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_consolidated_subchrom"
        ANNO_MERGED_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_subchrom"
        ANNOTATED="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_consolidated_${src_cell}"
        ANNO_MERGED="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_${src_cell}"
        if [[ ! -e $ANNO_MERGED ]]; then
            mkdir $ANNO_MERGED
        fi
        ANNOTATED_WITH_TSS_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_combined_subchrom"
        ANNOTATED_WITH_TSS="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_combined_${src_cell}"
        if [[ ! -e $ANNOTATED_WITH_TSS ]]; then
            mkdir $ANNOTATED_WITH_TSS
        fi  
        ANNOTATED_TSS_AND_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_${src_cell}_sub_consolidated_tss_and"
        ANNOTATED_TSS_AND="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_${src_cell}_consolidated_tss_and"
        if [[ ! -e $ANNOTATED_TSS_AND ]]; then
            mkdir $ANNOTATED_TSS_AND
        fi         
        ANNOTATED_PERM="$BASE/annotated_consolidated_perm_${src_cell}"
        ANNO_MERGED_PERM_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_perm_subchrom"
        ANNO_MERGED_PERM="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_perm_${src_cell}"
        if [[ ! -e $ANNO_MERGED_PERM ]]; then
            mkdir $ANNO_MERGED_PERM
        fi
        ANNO_MERGED_RPKM="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_rpkm_${src_cell}"
        if [[ ! -e $ANNO_MERGED_RPKM ]]; then
            mkdir $ANNO_MERGED_RPKM
        fi
        PRECISION_RECALL_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/precision_recall_subchrom"
        PRECISION_RECALL="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/precision_recall_${src_cell}"
        if [[ ! -e $PRECISION_RECALL ]]; then
            mkdir $PRECISION_RECALL
        fi
        UNKNOWNS_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/unknown_subchrom"
        UNKNOWNS="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/unknown_${src_cell}"
        if [[ ! -e $UNKNOWNS ]]; then
            mkdir $UNKNOWNS
        fi
        WIG="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/wig_chroms/${dest_cell}.chr"
        MISCLASSIFIED_PLOTS_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/misclassified_subchrom"
        MISCLASSIFIED_PLOTS="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/misclassified_${src_cell}"
        if [[ ! -e $MISCLASSIFIED_PLOTS ]]; then
            mkdir $MISCLASSIFIED_PLOTS
        fi
        TSS_MERGED_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/tss_merged_subchrom"
        TSS_MERGED="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/tss_merged_${src_cell}"
        if [[ ! -e $TSS_MERGED ]]; then
            mkdir $TSS_MERGED
        fi
        ANNOTATED_RPKM="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_rpkm"
        if [[ ! -e $ANNOTATED_RPKM ]]; then
            mkdir $ANNOTATED_RPKM
        fi
        ANNOTATED_RPKM_50BP="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_rpkm_50bp"
        if [[ ! -e $ANNOTATED_RPKM_50BP ]]; then
            mkdir $ANNOTATED_RPKM_50BP
        fi 
        ANNO_MERGED_RPKM="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_rpkm"
        if [[ ! -e $ANNO_MERGED_RPKM ]]; then
            mkdir $ANNO_MERGED_RPKM
        fi
        ANNOTATION_DISTRIBS_SUB="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotation_distributions_subchrom"
        ANNOTATION_DISTRIBS="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotation_distributions"
        if [[ ! -e $ANNOTATION_DISTRIBS ]]; then
            mkdir $ANNOTATION_DISTRIBS
        fi
        
        #Split the ChromHMM annotation file by chromosome. This only needs to be done once for each cell type.
        CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${dest_cell}_chromhmm_15_liftOver.bed"
        CHROMHMM_SPLIT="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${dest_cell}_split"
        if [[ ! -e $CHROMHMM_SPLIT ]]; then
            mkdir $CHROMHMM_SPLIT
            for chr in `cut -f 1 $CHROMHMM| sort | uniq`; do
                grep -w $chr $CHROMHMM > $CHROMHMM_SPLIT/$chr.bed
            done
        fi
        
        # Get the precision and recall for each chromosome, for:
        # 1. Our predictions
        # 2. TSS-based predictions
        # 3. Our predictions combined with TSS
        # 4. Permuted predictions
        # 5. Predictions using RPKM intensity
        for chr in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y';
            do       
                #Intersect TSS-based annotations with the ground truth.
                ./run_get_regions $WIG_CHROMS/${dest_cell}.chr${chr}.wig $TSS_PROMOTER_SRC/chr${chr}.bed ${chr} $TSS_PROMOTER_SHAPES/shapes${chr} $TSS_PROMOTER/final${chr}.bed
                bedtools intersect -wao -a $TSS_PROMOTER/final${chr}.bed -b $CHROMHMM > $TSS_MERGED/anno${chr}.bed
                
                #Intersect our annotations with the ground truth.
                bedtools intersect -wao -a $ANNOTATED/anno${chr}.bed -b $CHROMHMM > $ANNO_MERGED/anno${chr}.bed

                #Combine shape-based and TSS-based annotations (OR) and intersect with the ground truth.
                python ../annotation_scripts/combine_prediction_beds.py $ANNOTATED/anno${chr}.bed $TSS_PROMOTER/final${chr}.bed $ANNOTATED_WITH_TSS/${chr}.bed
                ./run_get_lines $WIG_CHROMS/${dest_cell}.chr${chr}.wig $ANNOTATED_WITH_TSS/${chr}.bed ${chr} $ANNOTATED_WITH_TSS/shapes${chr} $ANNOTATED_WITH_TSS/final${chr}.bed
                bedtools intersect -wao -a $ANNOTATED_WITH_TSS/${chr}.bed -b $CHROMHMM > $ANNO_MERGED/anno${chr}_tss.bed
                
                #Intersect shape AND TSS-based annotations with ground truth.
                bedtools intersect -wao -a $ANNOTATED_TSS_AND/anno${chr}.bed -b $CHROMHMM > $ANNO_MERGED/anno${chr}_tss_and.bed
                
                #Make RPKM-based predictions and intersect the RPKM-predicted data with the ground truth.
                #This only needs to be done once per cell line.
                if [[ ! -e $ANNOTATED_RPKM/anno${chr}.bed ]]; then
                    python ../annotation_scripts/predict_from_rpkm.py $WIG_CHROMS/${dest_cell}.chr${chr}.wig $ANNOTATED_RPKM_50BP/anno${chr}.bed $ANNOTATED_RPKM/anno${chr}.bed $chr
                    Intersect the RPKM data with the ground truth and retrieve the signal.
                    bedtools intersect -wao -a $ANNOTATED_RPKM/anno${chr}.bed -b $CHROMHMM > $ANNO_MERGED_RPKM/anno${chr}.bed
                    ./run_get_lines $WIG_CHROMS/${dest_cell}.chr${chr}.wig $ANNOTATED_RPKM/anno${chr}.bed ${chr} $ANNOTATED_RPKM/shapes_anno${chr} $ANNOTATED_RPKM/anno${chr}final.bed
                    awk {'printf ("%s\t%s\t%s\t%s\t1.0\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9)'} $ANNO_MERGED_RPKM/anno${chr}.bed > $ANNO_MERGED_RPKM/anno${chr}_final.bed
                fi 
                
                #Intersect the permuted data with the ground truth.
                bedtools intersect -wao -a $ANNOTATED_PERM/anno${chr}.bed -b $CHROMHMM > "$ANNO_MERGED_PERM/anno${chr}.bed"
            done
        
        #Plot precision and recall.    
        python plot_precision_recall.py $ANNO_MERGED/ $TSS_MERGED/ $ANNO_MERGED_RPKM/ $ANNO_MERGED_PERM_SUB/ $ANNOTATED/ $TSS_PROMOTER_SHAPES/ $ANNOTATED_WITH_TSS/ $ANNOTATED_TSS_AND/ $ANNOTATED_RPKM/ $ANNOTATED_PERM/ $BASE/  $PRECISION_RECALL/ $UNKNOWNS/ $WIG ${dest_cell} ${src_cell} 0
        
        if [ $i -eq 4 ]; then
            python plot_precision_recall.py $ANNO_MERGED_SUB/ $TSS_MERGED_SUB/ $ANNO_MERGED_RPKM/ $ANNO_MERGED_PERM/ $ANNOTATED_SUB/ $TSS_PROMOTER_SHAPES/ $ANNOTATED_WITH_TSS_SUB/ $ANNOTATED_TSS_AND_SUB/ $ANNOTATED_RPKM/ $ANNOTATED_PERM/ $BASE/sub  $PRECISION_RECALL_SUB/ $UNKNOWNS_SUB/ $WIG ${dest_cell} ${src_cell} 1
        fi
        
        #Plot the distribution of annotations for all methods and the ground truth.
        python plot_annotation_distribs.py $ANNOTATED $ANNOTATED_WITH_TSS $ANNOTATED_TSS_AND $CHROMHMM_SPLIT $ANNOTATED_PERM ${src_cell} ${dest_cell} $ANNOTATION_DISTRIBS/${src_cell}_to_${dest_cell}
                   
        #Plot the distribution of unknown regions
        python line_plot_unknown_distribs.py $ANNO_MERGED/anno $WIG $ANNOTATED/shapes_anno $UNKNOWNS $MISCLASSIFIED_PLOTS ${src_cell} ${dest_cell}
    done 