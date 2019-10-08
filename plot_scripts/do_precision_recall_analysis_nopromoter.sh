#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

#Directories that do not need to be created for each source-dest combination.
 SRC=(A549 A549 A549 Brain Brain Brain H1 H1 H1 GM12878)
 DEST=(A549 Brain H1 A549 Brain H1 A549 Brain H1 GM12878)

#Run analysis for each source and destination pair.
#for i in 0 1 2 3 4 5 6 7 8;
for i in 9;
    do
        src_cell=${SRC[$i]}
        dest_cell=${DEST[$i]}
        #Directories used in analysis
        BASE="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}"
        WIG="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/wig_chroms/${dest_cell}.chr"
        ANNOTATED="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_consolidated_nopromoter_05"
        ANNO_MERGED="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/annotated_merged_nopromoter_05"
        if [[ ! -e $ANNO_MERGED ]]; then
            mkdir $ANNO_MERGED
        fi
        PRECISION_RECALL="/fs/project/PAS0272/Tara/DNase_SOM/${dest_cell}/precision_recall_${src_cell}_nopromoter_05"
        if [[ ! -e $PRECISION_RECALL ]]; then
            mkdir $PRECISION_RECALL
        fi
        
        #CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/${dest_cell}_chromhmm_nopromoter.bed"
        CHROMHMM="/fs/project/PAS0272/Tara/DNase_SOM/chromHmm/GM_colored.bed_nopromoter.bed"
        
        # Get the precision and recall for each chromosome
        for chr in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y';
            do       
                
                #Intersect our annotations with the ground truth.
                bedtools intersect -wao -a $ANNOTATED/anno${chr}.bed -b $CHROMHMM > $ANNO_MERGED/anno${chr}.bed
                
            done
        
        #Plot precision and recall.
        for thresh in 1 2 5 10 20 30 40;
            do        
                python plot_precision_recall_nopromoter_abovethreshonly.py $ANNO_MERGED/ $ANNOTATED/ $BASE/  $PRECISION_RECALL/ $WIG ${dest_cell} ${src_cell} 0 $thresh
            done
        
    done 