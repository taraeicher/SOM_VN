#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#!/bin/bash   
#Need to install seaborn

SRC=""
BASE_FILENAME=""
ANNOTATED_ALL_REG=""
SHAPES_COMPREHENSIVE=""
SHAPES=""
LOG=""
while getopts n:d:s:m:l:r: option; do
        case "${option}" in
            n) SRC=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            s) SHAPES_COMPREHENSIVE=$(realpath $OPTARG);;
            m) SHAPES=$(realpath $OPTARG);;
            l) LOG=$(realpath $OPTARG);;
            r) ANNOTATED_ALL_REG=$(realpath $OPTARG);;
        esac
    done
enhancer=0
SOM="$BASE_FILENAME/$SRC/som_output_final"
VAL="$BASE_FILENAME/$SRC/cluster_validitytest" 
SHAPEFILE_FIGS="$BASE_FILENAME/$SRC/annotation_figs/${SRC}test"

#Create all needed directories.
if [[ ! -e $SHAPEFILE_FIGS ]]; then
	mkdir $SHAPEFILE_FIGS
fi

# Plot the cluster validity heatmap.
python compute_validity.py $ANNOTATED_ALL_REG/ $SOM/ $ANNOTATED_ALL_REG/ $VAL

python annotation_similarity_heatmap.py ${SHAPES_COMPREHENSIVE} ${LOG} ${SHAPEFILE_FIGS}_ratio ${SHAPEFILE_FIGS}_heatmap ${SHAPES} $enhancer $SRC
