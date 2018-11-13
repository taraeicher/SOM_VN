CELL="H1"
BASE="/fs/project/PAS0272/Tara/DNase_SOM"

#Get blacklist regions.
cd $BASE
#wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
#gunzip hg38.blacklist.bed.gz

#Separate blacklist into chromosomes.
mkdir -p blacklist_chroms_hg38
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; 
    do
        #Get the blacklist regions associated with the chromosome.
        #grep -w chr$chr hg38.blacklist.bed > blacklist_chroms_hg38/chr${chr}.bed
        
        #Get all 50 bp regions with RPKM greater than 50.
        awk ' $2 >= 50.0 ' ${BASE}/${CELL}/wig_chroms/${CELL}.chr${chr}.wig > ${CELL}.chr${chr}_over50.wig
        
        #Print percentage of 50 bp regions with RPKM greater than 50.
        over50=$(wc -l < "${CELL}.chr${chr}_over50.wig")
        all=$(wc -l < "${BASE}/${CELL}/wig_chroms/${CELL}.chr${chr}.wig")
        all=$(($all - 2))
        percentile=$(bc -l <<< "($all - $over50) / $all")
        echo $percentile
        
        #Overlap 50 bp regions with RPKM greater than 50 with blacklists. See how many overlap.
        count_blacklist=0
        for p in $(awk '{print $1}' Brain.chr${chr}_over50.wig);
            do
                if [ -f "blacklist_chroms_hg38/chr${chr}.bed" ]; then
                    new_count=$(awk '{if ($2 <= $p && $3 >= $p) print $2}' blacklist_chroms_hg38/chr${chr}.bed | wc -l)
                    count_blacklist=$(($count_blacklist + $new_count))
                fi
            done
            
        #Get the percentage of bins over 50 that overlap with blacklists.
        percent=0
        if [ $count_blacklist -gt 0 ]; then
            percent=$(bc -l <<< "$count_blacklist / $over50")
        fi
        echo $percent
    done
    