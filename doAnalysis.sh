#PBS -l nodes=2:ppn=28
#PBS -l walltime=48:00:00 
#!/bin/bash   
"""
Dependencies:
1. Following python modules for gosr: argparse, logging, sys, itertools, numpy, pysam
2. Tensorflow, which can be installed here: https://www.tensorflow.org/install/install_linux
3. Samtools: https://sourceforge.net/projects/samtools/files/
4. Bowtie2: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/
5. Picard: https://broadinstitute.github.io/picard/
6. Wig-split: https://github.com/taoliu/taolib/blob/master/Scripts/wig_split.py
7. Taolib directory: https://github.com/taoliu/taolib
6. Added to your path: gosr/bin, bowtie2, samtools, Picard, taolib, wig_split
7. Added to your python path: pythonmodules/gosr/lib
8. Make sure you are not using python3, which causes an error with the print statements in gosr and in wig_split
9. Make sure your version of python is compiled to use UCS2
10. gap statistic: https://github.com/minddrummer/gap

Usage:
You can start with either a FASTQ file and reference genome or a BAM file. This script will convert it into a single bp resolution
WIG file, then detect latent classes according to signal patterns in the WIG.

Optionally, you can overlay the latent classes with reference classes. To do this, you will need to supply a second BAM file and a
corresponding BED file with the annotations.

To modify the number of window sizes used, you will need to modify the C code and recompile as well as modifying the TensorFlow code.

It is necessary to run the following commands before running this script:
module load python/2.7
module load cuda/8.0.44
"""
#Enable job control.
	set -m

#Variables
	CELL_LINE="/fs/project/PAS0272/Tara/DNase_SOM"
	REF_GENOME="$HOME/hg38/hg38"
	BASE_FILENAME="A549"
	PICARD_PATH="$HOME/picard.jar"
	CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
	SAMTOOLS_FILTER=1540
	WIG_SPLIT_PATH="$HOME/taolib/Scripts/"
	BIN_SIZE=50
	PYTHON_VERSION=2.7
	CUDA_VERSION=8.0.44
	WINDOW_SIZE_COUNT=7
	FASTQ_DIR="/fs/project/PAS0272/Tara/DNase_SOM/A549"
	
	
#Create all needed directories.
	if [[ ! -e $BASE_PATH/$CELL_LINE/som_output ]]; then
		mkdir $BASE_PATH/$CELL_LINE/som_output
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/dataFiles ]]; then
		mkdir $BASE_PATH/$CELL_LINE/dataFiles
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/dataFiles_anno ]]; then
		mkdir $BASE_PATH/$CELL_LINE/dataFiles_anno
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/dataFiles_shifted ]]; then
		mkdir $BASE_PATH/$CELL_LINE/dataFiles_shifted
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/som_output_filtered ]]; then
		mkdir $BASE_PATH/$CELL_LINE/som_output_filtered
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/som_output_shifted ]]; then
		mkdir $BASE_PATH/$CELL_LINE/som_output_shifted
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/som_output_final ]]; then
		mkdir $BASE_PATH/$CELL_LINE/som_output_final
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/wig_chroms ]]; then
		mkdir $BASE_PATH/$CELL_LINE/wig_chroms
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/anno_beds ]]; then
		mkdir $BASE_PATH/$CELL_LINE/anno_beds
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/anno_beds_consolidated ]]; then
		mkdir $BASE_PATH/$CELL_LINE/anno_beds_consolidated
	fi
	if [[ ! -e $BASE_PATH/$CELL_LINE/anno_beds_final ]]; then
		mkdir $BASE_PATH/$CELL_LINE/anno_beds_final
	fi
	
#FASTQ alignment and processing.
"""
	echo -e "---------------------------------------Aligning to reference genome at $REF_GENOME.-----------------------------------"
    bowtie2 -x $REF_GENOME -U $FASTQ_DIR/$CELL_LINE.fastq -k 1  -S $BASE_PATH/${BASE_FILENAME}_trimmed.bam --threads 16 --trim3 16 
	echo -e "-------------------------------------------------Alignment complete.--------------------------------------------------\n"
	
	echo -e "\n--------------------------------------------------Sorting BAM file.---------------------------------------------------"
    samtools sort -o $BASE_PATH/${BASE_FILENAME}_sorted.bam -T $BASE_PATH/${BASE_FILENAME}_tmp $BASE_PATH/${BASE_FILENAME}_trimmed.bam
	echo -e "----------------------------------------------------Sorting complete.-------------------------------------------------\n"
	
	echo -e "\n------------------------------------------------Indexing sorted BAM file.---------------------------------------------"
    samtools index $BASE_PATH/${BASE_FILENAME}_sorted.bam
	echo -e "-------------------------------------------------Indexing complete.------------------------------------------------------\n"
	
	echo -e "\n--------------------------------------------Removing poor-quality reads.-----------------------------------------------"
    samtools view -bq 1 $BASE_PATH/${BASE_FILENAME}_sorted.bam > $BASE_PATH/${BASE_FILENAME}_mapq1.bam   
	samtools index $BASE_PATH/${BASE_FILENAME}_mapq1.bam
	echo -e "--------------------------------------------------Removal complete.----------------------------------------------------\n"
	
	echo -e "\n-----------------------------------Removing mitochondrial and unknown reads.-------------------------------------------"
	samtools view -b $BASE_PATH/${BASE_FILENAME}_mapq1.bam $CHROMS > $BASE_PATH/${BASE_FILENAME}_nochrM.bam
    samtools index $BASE_PATH/${BASE_FILENAME}_nochrM.bam
	echo -e "--------------------------------------------------Removal complete.-----------------------------------------------------\n"
	
	echo -e "\n----------------------------------------------Removing duplicate reads.---------------------------------------------------"
    java -jar $PICARD_PATH MarkDuplicates INPUT=$BASE_PATH/${BASE_FILENAME}_nochrM.bam OUTPUT=$BASE_PATH/${BASE_FILENAME}_marked.bam METRICS_FILE=$BASE_PATH/${BASE_FILENAME}_metrics.txt
    samtools view -F $SAMTOOLS_FILTER -b $BASE_PATH/${BASE_FILENAME}_marked.bam > $BASE_PATH/${BASE_FILENAME}_unsorted.bam
    samtools index $BASE_PATH/${BASE_FILENAME}_unsorted.bam

	samtools sort -o $BASE_PATH/${BASE_FILENAME}.bam -T $BASE_PATH/${BASE_FILENAME}_tmp $BASE_PATH/${BASE_FILENAME}_unsorted.bam 
	echo -e "-------------------------------------------------Duplicate reads removed.-------------------------------------------------\n"
"""	
#Output RPKM intensities in a WIG file.
	echo -e "\n----------------------------------------------------Creating WIG file.--------------------------------------------------"
	module load python/$PYTHON_VERSION
	module load cuda/$CUDA_VERSION
	gosr binbam -f 0 -n 1000 -t $CELL_LINE $BASE_PATH/$CELL_LINE.bam $BIN_SIZE $CELL_LINE > $BASE_PATH/$CELL_LINE/$CELL_LINE.wig
	echo -e "------------------------------------------------------WIG file complete.-------------------------------------------------\n"
	
#Split WIG files by chromosome.
	echo -e "\n---------------------------------------------Splitting WIG file by chromosome.-------------------------------------------"
	python $WIG_SPLIT_PATH/wig_split.py $BASE_PATH/$CELL_LINE/$CELL_LINE.wig $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE
	if [[ -e $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*chrM* ]]; then
		rm $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*chrM*
	fi
	if [[ -e $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*random* ]]; then
		rm $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*random*
	fi
	if [[ -e $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*chrM* ]]; then
		rm $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*chrM*
	fi
	if [[ -e $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*_g* ]]; then
		rm $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE/*_g*
	fi
	echo -e "-----------------------------------------------------Splitting complete.------------------------------------------------\n"
	
#Data preprocessing
	echo -e "\n---------------------------------------------Compiling data processing code.--------------------------------------------"
	gcc -pthread -lm -o runGetData getFileData.c
	echo -e "-----------------------------------------------------Compiling complete.------------------------------------------------\n"
	
#Method for running the pipeline for a chromosome.
	run_pipeline() {
		local file=$1
		c=${file#*wig_chroms/chr} 

		#Generate input files for training and annotation.
		echo -e "\n------------------------------Formatting data for training and annotation for chrom $c.-------------------------------------"
		./runGetData $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE.chr$c.wig $BIN_SIZE 0 Y $c $BASE_PATH/$CELL_LINE/dataFiles/chrom$c
		./runGetData $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE.chr$c.wig $BIN_SIZE 0 N $c $BASE_PATH/$CELL_LINE/dataFiles_anno/chrom$c
		echo -e "-------------------------------------Data formatting complete for chrom $c.------------------------------------------------\n"
		
		#Shuffle the input files.
		echo -e "\n--------------------------------------Shuffling the training data for chrom $c.--------------------------------------------"
		cd $BASE_PATH/$CELL_LINE/dataFiles/
		window="window"
		for i in $(seq 1 $WINDOW_SIZE_COUNT)
			do
				shuf chrom$c$window$i > shuf_chrom$c$window$i
			done
		echo -e "--------------------------------------------Shuffling complete for chrom $c.-------------------------------------------------\n"
		cd $BASE_PATH/$CELL_LINE
			
		#Shift the input to its best representation.
		echo -e "\n-------------------------------------Shifting the data to center peaks for chrom $c.-----------------------------------------"
		python shift_input.py $BASE_PATH/$CELL_LINE/dataFiles/shuf_chrom$c$window $BASE_PATH/$CELL_LINE/dataFiles_shifted/chrom$c$window $BIN_SIZE $BASE_PATH/$CELL_LINE/wig_chroms/$CELL_LINE.chr$c.wig
		echo -e "----------------------------------------------Shifting complete for chrom $c.-------------------------------------------------\n"
		
		#Run the SOM.
		echo -e "\n--------------------------------Building the SOM model for chrom $c. This may take a while.-----------------------------------"
		module load cuda/$CUDA_VERSION
		python som_auto.py $BASE_PATH/$CELL_LINE/dataFiles_shifted/chrom$c $BASE_PATH/$CELL_LINE/som_output/chrom$c $BASE_PATH/$CELL_LINE/$CELL_LINE.wig
		echo -e "---------------------------------------------SOM model is ready for chrom $c.---------------------------------------------------\n"
		
		#Remove all SOM centroids to which no regions map.
		echo -e "\n----------------------------------------Removing unused SOM centroids $c.--------------------------------------------"
		python remove_by_cutoff.py $BASE_PATH/$CELL_LINE/som_output/som_centroid_chrom$c $WINDOW_SIZE_COUNT 1 $BASE_PATH/$CELL_LINE/som_output_filtered/som_centroid_chrom$c
		echo -e "------------------------------------------------Removal complete for chrom $c.-----------------------------------------------------\n"
		
		#Merge shifted regions.
		echo -e "\n----------------------------------------Merging shifted centroids for chrom $c.--------------------------------------------------"
		python merge_shifted.py $WINDOW_SIZE_COUNT $BASE_PATH/$CELL_LINE/som_output_filtered/som_centroid_chrom$c $BASE_PATH/$CELL_LINE/som_output_shifted/som_centroid_chrom$c
		echo -e "------------------------------------------------Merging complete for chrom $c.------------------------------------------------------\n"
		
		#Remove duplicate centroids using kmeans.
		echo -e "\n--------------------------------Obtaining the final clusters using K-means for chrom $c.-----------------------------------------"
		python kmeans_centroids.py $WINDOW_SIZE_COUNT $BASE_PATH/$CELL_LINE/som_output_shifted/som_centroid_chrom$c $BASE_PATH/$CELL_LINE/som_output_final/som_centroid_chrom$c
		echo -e "-------------------------------------------------K-means complete for chrom $c.---------------------------------------------------\n"
		
		#Annotate regions with cluster data.
		echo -e "\n------------------------------------Annotating the regions by cluster for chrom $c.----------------------------------------------"
		python make_cluster_bed.py $BASE_PATH/$CELL_LINE/dataFiles_anno/chrom$c$window $BASE_PATH/$CELL_LINE/som_output_final/som_centroid_chrom$c $BASE_PATH/$CELL_LINE/anno_beds/anno$c
		echo -e "------------------------------------------------Annotation complete for chrom $c.------------------------------------------------\n"
	
		echo -e "\n------------------------------------Consolidating annotations per window size for chrom $c.-------------------------------------\n"
		if test "$(ls -A "$BASE_PATH/anno_beds/")"; then
			for bedfile in $BASE_PATH/$CELL_LINE/anno_beds/*;
				do
					pos=$(echo $bedfile | grep -b -o anno_beds/ | awk 'BEGIN {FS=":"}{print $1}')
					after=${bedfile:pos}
					pos=$pos+10
					after=$(echo ${bedfile:$pos})
					bedtools sort -i  $BASE_PATH/$CELL_LINE/anno_beds/$after > $BASE_PATH/$CELL_LINE/anno_beds/sorted_$after
				done
		fi
		python consolidate_each_window.py $BASE_PATH/$CELL_LINE/anno_beds/sorted_anno$c $BASE_PATH/$CELL_LINE/anno_beds_consolidated/anno$c
		echo -e "\n------------------------------------Consolidating per window size complete for chrom $c.-----------------------------------------"
		#python consolidate_all_windows.py $BASE_PATHanno_beds_consolidated/anno$c $BASE_PATHanno_beds_final/anno$c &
	}
	#Run the pipeline for each chromosome separately.
	if test "$(ls -A "$BASE_PATH/wig_chroms/")"; then
		for f in $BASE_PATH/$CELL_LINE/wig_chroms/*;
			do 
				run_pipeline $f & 
			done
	fi
	wait
	echo Done!
	exit 0
		#H3K27ac A549 https://www.encodeproject.org/experiments/ENCSR783SNV/
		#H3K27me3 H1 https://www.encodeproject.org/experiments/ENCSR186OBR/ (ONLY hg38)