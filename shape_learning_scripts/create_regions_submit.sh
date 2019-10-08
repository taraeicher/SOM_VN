#PBS -l nodes=1:ppn=28
#PBS -l walltime=2:00:00 
#!/bin/bash   
cd /fs/project/PAS0272/Tara/DNase_SOM/shape_learning_scripts/
./create_regions.sh -n Brain -d /fs/project/PAS0272/Tara/DNase_SOM/ -b /fs/project/PAS0272/Tara/DNase_SOM/bam/ -s=/fs/project/PAS0272/Tara/DNase_SOM/hg38.chrom.sizes -l /fs/project/PAS0272/Tara/DNase_SOM/hg38.blacklist_merged.bed