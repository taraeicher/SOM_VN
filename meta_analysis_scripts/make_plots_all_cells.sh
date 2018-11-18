#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

WIG_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/Brain/wig_chroms/"
WIG_A549="/fs/project/PAS0272/Tara/DNase_SOM/A549/wig_chroms/"
WIG_H1="/fs/project/PAS0272/Tara/DNase_SOM/H1/wig_chroms/"
TSS="/fs/project/PAS0272/Tara/DNase_SOM/TSS_hg38_dedupe.bed"
WIG_DISTRIB="/fs/project/PAS0272/Tara/DNase_SOM/wig_distribs"
CHROMHMM_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/Brain_chromhmm_15_liftOver.bed"
CHROMHMM_A549="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/A549_chromhmm_15_liftOver.bed"
CHROMHMM_H1="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/H1_chromhmm_15_liftOver.bed"
CHROMHMM_DISTRIB="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm_distribs"
BASE="/fs/project/PAS0272/Tara/DNase_SOM/"

#Plot all true distributions for each predicted distribution.
python plot_true_distribs_all.py $BASE $BASE/truedistribs

#Plot the distribution of RPKM signal per cell type.
#python plot_wig_distribs_violin.py $WIG_A549 $WIG_H1 $WIG_BRAIN $WIG_DISTRIB

#Plot the distribution of ChromHMM region widths.
#python plot_chromhmm_distribs_violin.py $CHROMHMM_A549 $CHROMHMM_H1 $CHROMHMM_BRAIN $CHROMHMM_DISTRIB

#Plot averaged precision and recall values in a single plot.
ANNO_MERGED="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_${src_cell}"
PRECISION_RECALL_START="/fs/project/PAS0272/Tara/DNase_SOM/"
PRECISION_RECALL_MIDDLE="/precision_recall_"
PRECISION_RECALL_END="/report_"
#python plot_precision_recall_all.py $PRECISION_RECALL_START $PRECISION_RECALL_MIDDLE $PRECISION_RECALL_END $BASE

#Plot averaged annotation distributions in a single plot.
ANNOTATION_DISTRIBS_START="/fs/project/PAS0272/Tara/DNase_SOM/"
ANNOTATION_DISTRIBS_END="/annotation_distributions"
#python plot_annotation_distribs_all.py $ANNOTATION_DISTRIBS_START $ANNOTATION_DISTRIBS_END ${BASE}/annotation_distribs
# #Print clusters annotated with annotation count and average TSS count.
for src in 'A549' 'Brain' 'H1';
    do  
        #Create directories.
        DATABASE="/fs/project/PAS0272/Tara/DNase_SOM/$src/database" 
        CLUSTER_PLOTS="/fs/project/PAS0272/Tara/DNase_SOM/$src/cluster_plots"
        TSS_COVERAGE_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/Brain/tss_coverage_${src}"
        TSS_COVERAGE_A549="/fs/project/PAS0272/Tara/DNase_SOM/A549/tss_coverage_${src}"
        TSS_COVERAGE_H1="/fs/project/PAS0272/Tara/DNase_SOM/H1/tss_coverage_${src}"
        ANNOTATED_BRAIN="/fs/project/PAS0272/Tara/DNase_SOM/Brain/annotated_consolidated_${src}"
        ANNOTATED_A549="/fs/project/PAS0272/Tara/DNase_SOM/A549/annotated_consolidated_${src}"
        ANNOTATED_H1="/fs/project/PAS0272/Tara/DNase_SOM/H1/annotated_consolidated_${src}"
        
        #Create all needed directories.
        if [[ ! -e $CLUSTER_PLOTS ]]; then
            mkdir $CLUSTER_PLOTS
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
        # for chr in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y'; 
            # do
                # bedtools coverage -counts -a $ANNOTATED_BRAIN/anno${chr}clust.bed -b $TSS > $TSS_COVERAGE_BRAIN/anno${chr}.bed
                # bedtools coverage -counts -a $ANNOTATED_A549/anno${chr}clust.bed -b $TSS > $TSS_COVERAGE_A549/anno${chr}.bed
                # bedtools coverage -counts -a $ANNOTATED_H1/anno${chr}clust.bed -b $TSS > $TSS_COVERAGE_H1/anno${chr}.bed
            # done
            
        #python print_annotated_clusters.py $DATABASE $CLUSTER_PLOTS/ $ANNOTATED_BRAIN/anno $ANNOTATED_A549/anno $ANNOTATED_H1/anno $TSS_COVERAGE_BRAIN/anno $TSS_COVERAGE_A549/anno $TSS_COVERAGE_H1/anno $src
    done 