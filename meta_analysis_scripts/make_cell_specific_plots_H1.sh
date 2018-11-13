#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

SRC="H1"
SOM="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/som_output_final"
#DATABASE_ALL="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/percentage"
#DATABASE="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/chromHmm_merged"
#LOG="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/ig_chromHmm_log"
DATABASE_ALL="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/database"
DATABASE="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/database_all"
LOG="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/database_log"
VAL="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/cluster_validity" 
SHAPEFILE_FIGS="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/annotation_figs/"
ANNOTATED_ALL_WIN="/fs/project/PAS0272/Tara/DNase_SOM/$SRC/anno_beds"

#Create all needed directories.
if [[ ! -e $SHAPEFILE_FIGS ]]; then
	mkdir $SHAPEFILE_FIGS
fi

# Plot the cluster validity heatmap.
#python compute_validity.py $ANNOTATED_ALL_WIN/ $SOM/ $ANNOTATED_ALL_WIN/ $VAL

#Plot the database similarity heatmap for each window size.
# for window in '2' '3' '4' '5' '6';
	# do 
        # python annotation_similarity_heatmap.py ${DATABASE_ALL}${window} ${LOG}${window} ${SHAPEFILE_FIGS}ratio${window} ${SHAPEFILE_FIGS}heatmap${window} ${window} ${DATABASE}${window}
	# done
cd "/fs/project/PAS0272/Tara/DNase_SOM/scripts"
python annotation_similarity_heatmap.py ${DATABASE_ALL} ${LOG} ${SHAPEFILE_FIGS} ${SHAPEFILE_FIGS} ${DATABASE}
