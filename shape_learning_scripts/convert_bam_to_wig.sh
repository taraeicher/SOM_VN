#PBS -l nodes=1:ppn=24
#PBS -l walltime=10:00:00 
#!/bin/bash   
# """
# Dependencies:
# 1. Tensorflow, which can be installed here: https://www.tensorflow.org/install/install_linux
# 2. Samtools: https://sourceforge.net/projects/samtools/files/
# 3. Taolib directory: https://github.com/taoliu/taolib
# 4. Gap statistic code: https://github.com/minddrummer/gap
# 5. Added to your path: gap
# 6. Bedtools
# 7. Kent utilities
# 8. Bamtools
# 9. BamCoverage

#Variables
    USAGE="This script is used for creating training regions from an input BAM file. It first converts the BAM file to a set of BAM files, one for each chromosome. Then, it converts these files to WIG files containing a single RPKM value at each 50 bp bin.\n
    <-n> The name of the cell line (e.g. Brain)\n
    <-d> The base filename where the input and output files will be stored (e.g. '/root/annoshaperun/').\n
    <-b> The BAM file used as input.\n
    <-i> The bin size used to generate the WIG file (default: 50 bp)\n
    <-s> The file containing a list of chromosome sizes. This is needed for splitting the BAM file by chromosome.\n
    <-l> The blacklist regions to exclude."
    
    CELL_LINE=""
    BASE_FILENAME=""
    BAM=""
    CHROMS_NUM="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    BIN_SIZE=50
    BLACKLIST=""
    while getopts n:d:b:i:s:l: option; do
        case "${option}" in
            n) CELL_LINE=$OPTARG;;
            d) BASE_FILENAME=$(realpath $OPTARG);;
            b) BAM=$(realpath $OPTARG);;
            i) BIN_SIZE=$OPTARG;;
            s) CHROMSIZES=$OPTARG;;
            l) BLACKLIST=$OPTARG;;
        esac
    done
    BASE_PATH=$BASE_FILENAME/$CELL_LINE
    
#Print message to user.
    echo -e $USAGE
    echo "Creating regions using the following settings:"
    echo "Name is $CELL_LINE"
    echo "Base directory is $BASE_FILENAME"
    echo "BAM file is $BAM"
    echo "Bin size in WIG file is $BIN_SIZE"
    echo -e "-------------------------------------------------------------------------------------------------------------------------\n"
	
#Create all needed directories.
    if [[ ! -e ${BASE_PATH}/bed_chroms ]]; then
        mkdir ${BASE_PATH}/bed_chroms
    fi
    if [[ ! -e ${BASE_PATH}/bigwig ]]; then
        mkdir ${BASE_PATH}/bigwig
    fi
    if [[ ! -e ${BASE_PATH}/bedgraph ]]; then
        mkdir ${BASE_PATH}/bedgraph
    fi
    if [[ ! -e ${BASE_PATH}/wig ]]; then
        mkdir ${BASE_PATH}/wig
    fi
    
    convert_to_wig() {
        chr=$1
        bedtools genomecov -ibam $BAM/$CELL_LINE.REF_chr$chr.bam -bg > ${BASE_PATH}/bedgraph/$chr.bg
        bamtools sort -in $BAM/$CELL_LINE.REF_chr$chr.bam -out $BAM/$CELL_LINE.REF_chr${chr}_sorted.bam 
        python create_index_pysam.py $BAM/$CELL_LINE.REF_chr${chr}_sorted.bam
        /users/PAS0272/osu5316/.local/bin/bamCoverage -b $BAM/$CELL_LINE.REF_chr${chr}_sorted.bam -o ${BASE_PATH}/bigwig/$chr.bw -bs $BIN_SIZE -bl $BLACKLIST --normalizeUsing RPKM
        bigWigToWig ${BASE_PATH}/bigwig/$chr.bw ${BASE_PATH}/wig/${chr}_unfiltered.wig
        awk -v chrom="chr${chr}" '{ if ($1 == chrom) { print } }' ${BASE_PATH}/wig/${chr}_unfiltered.wig > ${BASE_PATH}/wig/${chr}.wig
        rm ${BASE_PATH}/wig/${chr}_unfiltered.wig
    }
    
    module load bamtools/2.2.2
    bamtools split -in $BAM/$CELL_LINE.bam -reference 
    pids=""
    for c in $CHROMS_NUM;
        do 
            convert_to_wig $c & 
            pids="$pids $!"
        done
     
    wait $pids
    echo "All conversions from BAM to WIG complete!"