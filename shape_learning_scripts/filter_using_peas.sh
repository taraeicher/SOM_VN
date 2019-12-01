#PBS -l nodes=1:ppn=2
#PBS -l walltime=2:00:00 

#Variables
CHROMHMM=""
CHROMS_NUM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
WIG=""
TRAINING_FILES=""
RESOLUTION=10
REGION_SZ=1000
while getopts b:c:r:t:w: option; do
   case "${option}" in
       b) RESOLUTION=$OPTARG;;
       c) CHROMHMM=$(realpath $OPTARG);;
       r) REGION_SZ=$OPTARG;;
       t) TRAINING_FILES=$(realpath $OPTARG);;
       w) WIG=$WIG;;
   esac
done
PYTHON_VERSION=3.6	
module load python/$PYTHON_VERSION

if [[ ! -e $TRAINING_FILES/peas ]]; then
    mkdir $TRAINING_FILES/peas
fi
if [[ ! -e $TRAINING_FILES/peas/percentile_cutoffs ]]; then
    mkdir $TRAINING_FILES/peas/percentile_cutoffs
fi
if [[ ! -e $WIG/peas ]]; then
    mkdir $WIG/peas
fi
	
# Method for running the pipeline for a chromosome.
run_pipeline() {
    local c=$1
    
    #Remove promoters from ChromHMM.
    awk '{ if ( $4 != "AP" && $4 != "OP" && $4 != "GE" && $4 != "TS") { print; } }' $CHROMHMM > ${CHROMHMM}_nopromoter.bed
    awk '{ if ( $4 == "AP" || $4 == "OP" || $4 == "GE" || $4 == "TS") { print; } }' $CHROMHMM > ${CHROMHMM}_onlypromoter.bed
    
    #Retain only training regions contained within an annotated ChromHMM region.
    python filter_regions.py $TRAINING_FILES/${c}.pkl ${CHROMHMM}_onlypromoter.bed $TRAINING_FILES/peas/${c}.pkl
    
    # Retain only WIG bins contained within an annotated ChromHMM region.
    python filter_wig.py $WIG/${c}.wig ${CHROMHMM}_onlypromoter.bed $WIG/peas/${c}.wig $c
    
    # Get the new intensity percentile of only the filtered wig.
    python ../common_scripts/get_intensity_percentile_wig.py $WIG/peas/${c}.wig $RESOLUTION $REGION_SZ 0.95 $TRAINING_FILES/peas/percentile_cutoffs/${c}.txt
}

#Run the pipeline from split WIG files to final set of annotations.
pids=""
for f in $CHROMS_NUM;
do 
    run_pipeline $f & 
    pids="$pids $!"
done
    