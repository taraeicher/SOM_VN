#K562 BAM - https://www.encodeproject.org/experiments/ENCSR921NMD/
wget https://www.encodeproject.org/files/ENCFF538GKX/@@download/ENCFF538GKX.bam
wget https://www.encodeproject.org/files/ENCFF299TOZ/@@download/ENCFF299TOZ.bam
cat ENCFF538GKX.bam ENCFF299TOZ.bam > K562.bam

#A549 BAM - https://www.encodeproject.org/experiments/ENCSR000ELW/
#NOTE: https://www.encodeproject.org/experiments/ENCSR136DNA/ did not work with gosr binbam.
#All positions were set to 0 in the WIG file.
wget https://www.encodeproject.org/files/ENCFF566SPI/@@download/ENCFF566SPI.bam
wget https://www.encodeproject.org/files/ENCFF414MBW/@@download/ENCFF414MBW.bam
cat ENCFF566SPI.bam ENCFF414MBW.bam > A549.bam
samtools sort -o A549_sorted.bam -T A549_tmp.bam A549.bam

#IMR-90 BAM - https://www.encodeproject.org/experiments/ENCSR477RTP/
wget https://www.encodeproject.org/files/ENCFF023JFG/@@download/ENCFF023JFG.bam
cp ENCFF023JFG.bam IMR90.bam

#H1-hESC BAM - https://www.encodeproject.org/experiments/ENCSR794OFW/
wget https://www.encodeproject.org/files/ENCFF059BEU/@@download/ENCFF059BEU.bam
cp ENCFF059BEU.bam H1.bam

#Brain tissue BAM - https://www.encodeproject.org/experiments/ENCSR595CSH/
wget https://www.encodeproject.org/files/ENCFF226FCY/@@download/ENCFF226FCY.bam
wget https://www.encodeproject.org/files/ENCFF982IRZ/@@download/ENCFF982IRZ.bam
cat ENCFF226FCY.bam ENCFF982IRZ.bam > Brain.bam
samtools sort -o Brain_sorted.bam -T Brain_tmp.bam Brain.bam

#K562 histone marks
#H3K4me3 - https://www.encodeproject.org/experiments/ENCSR668LDD/
wget https://www.encodeproject.org/files/ENCFF999RPS/@@download/ENCFF999RPS.bed.gz
gunzip ENCFF999RPS.bed.gz
cp ENCFF999RPS.bed K562_h3k4me3.bed
bedtools sort -i K562_h3k4me3.bed > K562_h3k4me3_sorted.bed 
#H3K27me3 - https://www.encodeproject.org/experiments/ENCSR000EWB/
wget https://www.encodeproject.org/files/ENCFF439XWU/@@download/ENCFF439XWU.bed.gz
gunzip ENCFF439XWU.bed.gz
cp ENCFF439XWU.bed K562_h3k27me3.bed
bedtools sort -i K562_h3k27me3.bed > K562_h3k27me3_sorted.bed 
#H3K27ac - https://www.encodeproject.org/experiments/ENCSR000AKP/
wget https://www.encodeproject.org/files/ENCFF437DPT/@@download/ENCFF437DPT.bed.gz
gunzip ENCFF437DPT.bed.gz
cp ENCFF437DPT.bed K562_h3k27ac.bed
bedtools sort -i K562_h3k27ac.bed > K562_h3k27ac_sorted.bed 

#A549 histone marks
#H3K27ac - https://www.encodeproject.org/experiments/ENCSR778NQS/
wget https://www.encodeproject.org/files/ENCFF932ORM/@@download/ENCFF932ORM.bed.gz
wget https://www.encodeproject.org/files/ENCFF637BNE/@@download/ENCFF637BNE.bed.gz
wget https://www.encodeproject.org/files/ENCFF697WAE/@@download/ENCFF697WAE.bed.gz
gunzip ENCFF932ORM.bed.gz ENCFF637BNE.bed.gz ENCFF697WAE.bed.gz
cat ENCFF932ORM.bed ENCFF637BNE.bed ENCFF697WAE.bed > A549_h3k27ac.bed
bedtools sort -i A549_h3k27ac.bed > A549_h3k27ac_sorted.bed 
bedtools merge -i A549_h3k27ac_sorted.bed > A549_h3k27ac_merged.bed
#H3K4me3 - https://www.encodeproject.org/experiments/ENCSR203XPU/
wget https://www.encodeproject.org/files/ENCFF404REU/@@download/ENCFF404REU.bed.gz
wget https://www.encodeproject.org/files/ENCFF956FPT/@@download/ENCFF956FPT.bed.gz
wget https://www.encodeproject.org/files/ENCFF857ZRY/@@download/ENCFF857ZRY.bed.gz
gunzip ENCFF404REU.bed.gz ENCFF956FPT.bed.gz ENCFF857ZRY.bed.gz
cat ENCFF404REU.bed ENCFF956FPT.bed ENCFF857ZRY.bed > A549_h3k4me3.bed
bedtools sort -i A549_h3k4me3.bed > A549_h3k4me3_sorted.bed 
bedtools merge -i A549_h3k4me3_sorted.bed > A549_h3k4me3_merged.bed

#IMR90 histone marks
#H3K27ac - https://www.encodeproject.org/experiments/ENCSR002YRE/
wget https://www.encodeproject.org/files/ENCFF732SKL/@@download/ENCFF732SKL.bed.gz
gunzip ENCFF732SKL.bed.gz 
cp ENCFF732SKL.bed IMR90_h3k27ac.bed
bedtools sort -i IMR90_h3k27ac.bed > IMR90_h3k27ac_sorted.bed 
#H3K27me3 - https://www.encodeproject.org/experiments/ENCSR431UUY/
wget https://www.encodeproject.org/files/ENCFF325CLJ/@@download/ENCFF325CLJ.bed.gz
gunzip ENCFF325CLJ.bed.gz
cp ENCFF325CLJ.bed IMR90_h3k27me3.bed
bedtools sort -i IMR90_h3k27me3.bed > IMR90_h3k27me3_sorted.bed 
#H3K4me3 - https://www.encodeproject.org/experiments/ENCSR087PFU/
wget https://www.encodeproject.org/files/ENCFF018CAH/@@download/ENCFF018CAH.bed.gz
gunzip ENCFF018CAH.bed.gz
cp ENCFF018CAH.bed IMR90_h3k4me3.bed
bedtools sort -i IMR90_h3k4me3.bed > IMR90_h3k4me3_sorted.bed 

#H1-hESC histone marks
#H3K27me3 - https://www.encodeproject.org/experiments/ENCSR928HYM/
wget https://www.encodeproject.org/files/ENCFF810VGL/@@download/ENCFF810VGL.bed.gz
gunzip ENCFF810VGL.bed.gz
cp ENCFF810VGL.bed H1_h3k27me3.bed
awk {'printf ("%s\t%s\t%s\thh3k27me3\n", $1, $2, $3)'} H1_h3k27me3.bed > H1_h3k27me3_piece.bed
#H3K4me3 - https://www.encodeproject.org/experiments/ENCSR443YAS/
wget https://www.encodeproject.org/files/ENCFF958TAD/@@download/ENCFF958TAD.bed.gz
gunzip ENCFF958TAD.bed.gz
cp ENCFF958TAD.bed H1_h3k4me3.bed
awk {'printf ("%s\t%s\t%s\thh3k4me3\n", $1, $2, $3)'} H1_h3k4me3.bed > H1_h3k4me3_piece.bed
#H3K27ac - https://www.encodeproject.org/experiments/ENCSR880SUY/
wget https://www.encodeproject.org/files/ENCFF714VTU/@@download/ENCFF714VTU.bed.gz
gunzip ENCFF714VTU.bed.gz
cp ENCFF714VTU.bed H1_h3k27ac.bed
awk {'printf ("%s\t%s\t%s\thh3k27ac\n", $1, $2, $3)'} H1_h3k27ac.bed > H1_h3k27ac_piece.bed
cat H1_h3k27me3_piece.bed H1_h3k4me3_piece.bed H1_h3k27ac_piece.bed > H1_histone_unsorted.bed
bedtools sort -i H1_histone_unsorted.bed > H1_histone.bed
bedtools merge -o distinct -c 4 -i H1_histone.bed > H1_histone_merged.bed

#Brain tissue histone marks
#H3K4me3 - https://www.encodeproject.org/experiments/ENCSR780FXX/
wget https://www.encodeproject.org/files/ENCFF297ZIH/@@download/ENCFF297ZIH.bed.gz
gunzip ENCFF297ZIH.bed.gz 
cp ENCFF297ZIH.bed Brain_h3k4me3.bed
bedtools sort -i Brain_h3k4me3.bed > Brain_h3k4me3_sorted.bed
awk {'printf ("%s\t%s\t%s\th3k4me3\n", $1, $2, $3)'} Brain_h3k4me3.bed > Brain_h3k4me3_piece.bed
#H3K27me3 - https://www.encodeproject.org/experiments/ENCSR997YTW/
wget https://www.encodeproject.org/files/ENCFF928EQD/@@download/ENCFF928EQD.bed.gz
gunzip ENCFF928EQD.bed.gz
cp ENCFF928EQD.bed Brain_h3k27me3.bed
bedtools sort -i Brain_h3k27me3.bed > Brain_h3k27me3_sorted.bed
awk {'printf ("%s\t%s\t%s\th3k27me3\n", $1, $2, $3)'} Brain_h3k27me3.bed > Brain_h3k27me3_piece.bed
cat Brain_h3k27me3_piece.bed Brain_h3k4me3_piece.bed > Brain_histone_unsorted.bed
bedtools sort -i Brain_histone_unsorted.bed > Brain_histone.bed
bedtools merge -o distinct -c 4 -i Brain_histone.bed > Brain_histone_merged.bed
