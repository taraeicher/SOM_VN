#GM12878 BAM - https://www.encodeproject.org/experiments/ENCSR000EMT/
wget https://www.encodeproject.org/files/ENCFF593WBR/@@download/ENCFF593WBR.bam
wget https://www.encodeproject.org/files/ENCFF445PZT/@@download/ENCFF445PZT.bam
cat ENCFF593WBR.bam ENCFF445PZT.bam > GM12878.bam
samtools sort -o GM12878_sorted.bam -T GM12878_tmp.bam GM12878.bam

#A549 BAM - https://www.encodeproject.org/experiments/ENCSR000ELW/
#NOTE: https://www.encodeproject.org/experiments/ENCSR136DNA/ did not work with gosr binbam.
#All positions were set to 0 in the WIG file.
wget https://www.encodeproject.org/files/ENCFF566SPI/@@download/ENCFF566SPI.bam
wget https://www.encodeproject.org/files/ENCFF414MBW/@@download/ENCFF414MBW.bam
cat ENCFF566SPI.bam ENCFF414MBW.bam > A549.bam
samtools sort -o A549_sorted.bam -T A549_tmp.bam A549.bam

#H1-hESC BAM - https://www.encodeproject.org/experiments/ENCSR794OFW/
wget https://www.encodeproject.org/files/ENCFF059BEU/@@download/ENCFF059BEU.bam
cp ENCFF059BEU.bam H1.bam

#Brain tissue BAM - https://www.encodeproject.org/experiments/ENCSR595CSH/
wget https://www.encodeproject.org/files/ENCFF226FCY/@@download/ENCFF226FCY.bam
wget https://www.encodeproject.org/files/ENCFF982IRZ/@@download/ENCFF982IRZ.bam
cat ENCFF226FCY.bam ENCFF982IRZ.bam > Brain.bam
samtools sort -o Brain_sorted.bam -T Brain_tmp.bam Brain.bam
