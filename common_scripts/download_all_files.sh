#GM12878 BAM
wget https://www.encodeproject.org/files/ENCFF593WBR/@@download/ENCFF593WBR.bam
wget https://www.encodeproject.org/files/ENCFF445PZT/@@download/ENCFF445PZT.bam
cat ENCFF593WBR.bam ENCFF445PZT.bam > GM12878.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E116_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E116_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E116_15_coreMarks_hg38lift_mnemonics.bed GM12878_E116.bed

#A549 BAM
wget https://www.encodeproject.org/files/ENCFF566SPI/@@download/ENCFF566SPI.bam
wget https://www.encodeproject.org/files/ENCFF414MBW/@@download/ENCFF414MBW.bam
cat ENCFF566SPI.bam ENCFF414MBW.bam > A549.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E114_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E114_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E114_15_coreMarks_hg38lift_mnemonics.bed A549_E114.bed

#H1-hESC BAM (low)
wget https://www.encodeproject.org/files/ENCFF059BEU/@@download/ENCFF059BEU.bam
cp ENCFF059BEU.bam H1.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E003_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E003_15_coreMarks_hg38lift_mnemonics.bed H1_E003.bed

#Brain tissue BAM (high)
wget https://www.encodeproject.org/files/ENCFF226FCY/@@download/ENCFF226FCY.bam
wget https://www.encodeproject.org/files/ENCFF982IRZ/@@download/ENCFF982IRZ.bam
cat ENCFF226FCY.bam ENCFF982IRZ.bam > Brain.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E081_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E081_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E081_15_coreMarks_hg38lift_mnemonics.bed Brain_E081.bed

#Brain tissue BAM (low)
wget https://www.encodeproject.org/files/ENCFF255UMI/@@download/ENCFF255UMI.bam
mv ENCFF255UMI.bam Brain_low.bam

#H1 BAM (high)
wget https://www.encodeproject.org/files/ENCFF915LEP/@@download/ENCFF915LEP.bam
mv ENCFF915LEP.bam H1_high.bam

#HeLa BAM (high)
wget https://www.encodeproject.org/files/ENCFF912JKA/@@download/ENCFF912JKA.bam
mv ENCFF912JKA.bam HeLa_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E117_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E117_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E117_15_coreMarks_hg38lift_mnemonics.bed HeLa_E117.bed

#HeLa BAM (low)
wget https://www.encodeproject.org/files/ENCFF342EWC/@@download/ENCFF342EWC.bam
mv ENCFF342EWC.bam HeLa_low.bam

#Heart tissue BAM (high)
wget https://www.encodeproject.org/files/ENCFF626HEO/@@download/ENCFF626HEO.bam
mv ENCFF626HEO.bam heart_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E083_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E083_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E083_15_coreMarks_hg38lift_mnemonics.bed heart_E083.bed

#Heart tissue BAM (low)
wget https://www.encodeproject.org/files/ENCFF837QAV/@@download/ENCFF837QAV.bam
mv ENCFF837QAV.bam heart_low.bam

#Stomach tissue BAM (high)
wget https://www.encodeproject.org/files/ENCFF696XWN/@@download/ENCFF696XWN.bam
mv ENCFF696XWN.bam stomach_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E092_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E092_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E092_15_coreMarks_hg38lift_mnemonics.bed stomach_E092.bed

# Stomach tissue BAM (low)
wget https://www.encodeproject.org/files/ENCFF899HTN/@@download/ENCFF899HTN.bam
mv ENCFF899HTN.bam stomach_low.bam

#B-Cell BAM (high)
wget https://www.encodeproject.org/files/ENCFF568CIR/@@download/ENCFF568CIR.bam
mv ENCFF568CIR.bam stomach_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E032_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E032_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E032_15_coreMarks_hg38lift_mnemonics.bed b_cell_E032.bed

#B-Cell BAM (low)
wget https://www.encodeproject.org/files/ENCFF799QKF/@@download/ENCFF799QKF.bam
mv ENCFF799QKF.bam b_cell_low.bam