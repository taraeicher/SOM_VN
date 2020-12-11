#GM12878 BAM
wget https://www.encodeproject.org/files/ENCFF593WBR/@@download/ENCFF593WBR.bam
wget https://www.encodeproject.org/files/ENCFF445PZT/@@download/ENCFF445PZT.bam
cat ENCFF593WBR.bam ENCFF445PZT.bam > GM12878.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E116_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E116_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E116_15_coreMarks_hg38lift_mnemonics.bed GM12878_E116.bed
#Peaks
wget https://www.encodeproject.org/files/ENCFF598KWZ/@@download/ENCFF598KWZ.bed.gz
gunzip ENCFF598KWZ.bed.gz
mv ENCFF598KWZ.bed GM12878.bed

#A549 BAM
wget https://www.encodeproject.org/files/ENCFF566SPI/@@download/ENCFF566SPI.bam
wget https://www.encodeproject.org/files/ENCFF414MBW/@@download/ENCFF414MBW.bam
cat ENCFF566SPI.bam ENCFF414MBW.bam > A549.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E114_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E114_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E114_15_coreMarks_hg38lift_mnemonics.bed A549_E114.bed
#Peaks
wget https://www.encodeproject.org/files/ENCFF529HMB/@@download/ENCFF529HMB.bed.gz
wget https://www.encodeproject.org/files/ENCFF823UOG/@@download/ENCFF823UOG.bed.gz
gunzip ENCFF529HMB.bed.gz ENCFF823UOG.bed.gz
cat ENCFF529HMB.bed ENCFF823UOG.bed > A549_cat.bed
sort -k1,1 -k2,2n A549_cat.bed > A549_cat_sorted.bed
bedtools merge -i A549_cat_sorted.bed > A549.bed

#H1-hESC BAM (low)
wget https://www.encodeproject.org/files/ENCFF059BEU/@@download/ENCFF059BEU.bam
cp ENCFF059BEU.bam H1.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E003_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E003_15_coreMarks_hg38lift_mnemonics.bed H1_E003.bed
#Peaks
wget https://www.encodeproject.org/files/ENCFF905XDS/@@download/ENCFF905XDS.bed.gz
gunzip ENCFF905XDS.bed.gz
mv ENCFF905XDS.bed H1_low.bed

#Brain tissue BAM (high)
wget https://www.encodeproject.org/files/ENCFF226FCY/@@download/ENCFF226FCY.bam
wget https://www.encodeproject.org/files/ENCFF982IRZ/@@download/ENCFF982IRZ.bam
cat ENCFF226FCY.bam ENCFF982IRZ.bam > Brain.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E081_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E081_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E081_15_coreMarks_hg38lift_mnemonics.bed Brain_E081.bed
#Peaks
wget https://www.encodeproject.org/files/ENCFF018ATG/@@download/ENCFF018ATG.bed.gz
wget https://www.encodeproject.org/files/ENCFF936VAD/@@download/ENCFF936VAD.bed.gz
gunzip ENCFF018ATG.bed.gz ENCFF936VAD.bed.gz
cat ENCFF018ATG.bed ENCFF936VAD.bed > Brain_high_cat.bed
sort -k1,1 -k2,2n Brain_high_cat.bed > Brain_high_cat_sorted.bed
bedtools merge -i Brain_high_cat_sorted.bed > Brain_high.bed

#Brain tissue BAM (low)
wget https://www.encodeproject.org/files/ENCFF255UMI/@@download/ENCFF255UMI.bam
mv ENCFF255UMI.bam Brain_low.bam
#Peaks
wget https://www.encodeproject.org/files/ENCFF149UIY/@@download/ENCFF149UIY.bed.gz
gunzip ENCFF149UIY.bed.gz
mv ENCFF149UIY.bed Brain_low.bed

#H1 BAM (high)
wget https://www.encodeproject.org/files/ENCFF915LEP/@@download/ENCFF915LEP.bam
mv ENCFF915LEP.bam H1_high.bam
#Peaks
wget https://www.encodeproject.org/files/ENCFF896QPX/@@download/ENCFF896QPX.bed.gz
gunzip ENCFF896QPX.bed.gz
mv ENCFF896QPX.bed H1_high.bed

#HeLa BAM (high)
wget https://www.encodeproject.org/files/ENCFF912JKA/@@download/ENCFF912JKA.bam
mv ENCFF912JKA.bam HeLa_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E117_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E117_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E117_15_coreMarks_hg38lift_mnemonics.bed HeLa_E117.bed
#Peaks
wget https://www.encodeproject.org/files/ENCFF950NDW/@@download/ENCFF950NDW.bed.gz
gunzip ENCFF950NDW.bed.gz
mv ENCFF950NDW.bed HeLa_high.bed

#HeLa BAM (low)
wget https://www.encodeproject.org/files/ENCFF342EWC/@@download/ENCFF342EWC.bam
mv ENCFF342EWC.bam HeLa_low.bam
#Peaks
wget https://www.encodeproject.org/files/ENCFF888QGO/@@download/ENCFF888QGO.bed.gz
gunzip ENCFF888QGO.bed.gz
mv ENCFF888QGO.bed HeLa_low.bed

#Heart tissue BAM (high)
wget https://www.encodeproject.org/files/ENCFF626HEO/@@download/ENCFF626HEO.bam
mv ENCFF626HEO.bam heart_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E083_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E083_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E083_15_coreMarks_hg38lift_mnemonics.bed heart_E083.bed
#Peaks
wget https://www.encodeproject.org/files/ENCFF911WSI/@@download/ENCFF911WSI.bed.gz
gunzip ENCFF911WSI.bed.gz
mv ENCFF911WSI.bed heart_high.bed

#Heart tissue BAM (low)
wget https://www.encodeproject.org/files/ENCFF837QAV/@@download/ENCFF837QAV.bam
mv ENCFF837QAV.bam heart_low.bam
#Peaks
wget https://www.encodeproject.org/files/ENCFF850WOE/@@download/ENCFF850WOE.bed.gz
gunzip ENCFF850WOE.bed.gz
mv ENCFF850WOE.bed heart_low.bed

#Stomach tissue BAM (high)
wget https://www.encodeproject.org/files/ENCFF696XWN/@@download/ENCFF696XWN.bam
mv ENCFF696XWN.bam stomach_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E092_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E092_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E092_15_coreMarks_hg38lift_mnemonics.bed stomach_E092.bed
# Peaks
wget https://www.encodeproject.org/files/ENCFF009YJE/@@download/ENCFF009YJE.bed.gz
gunzip ENCFF009YJE.bed.gz
mv ENCFF009YJE.bed stomach_high.bed

# Stomach tissue BAM (low)
wget https://www.encodeproject.org/files/ENCFF899HTN/@@download/ENCFF899HTN.bam
mv ENCFF899HTN.bam stomach_low.bam
# Peaks
wget https://www.encodeproject.org/files/ENCFF751PUA/@@download/ENCFF751PUA.bed.gz
gunzip ENCFF751PUA.bed.gz
mv ENCFF751PUA.bed stomach_low.bed

#B-Cell BAM (high)
wget https://www.encodeproject.org/files/ENCFF568CIR/@@download/ENCFF568CIR.bam
mv ENCFF568CIR.bam stomach_high.bam
#ChromHMM
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E032_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E032_15_coreMarks_hg38lift_mnemonics.bed.gz
mv E032_15_coreMarks_hg38lift_mnemonics.bed b_cell_E032.bed
# Peaks
wget https://www.encodeproject.org/files/ENCFF009YJE/@@download/ENCFF009YJE.bed.gz
gunzip ENCFF009YJE.bed.gz
mv ENCFF009YJE.bed b_cell_high.bed

#B-Cell BAM (low)
wget https://www.encodeproject.org/files/ENCFF799QKF/@@download/ENCFF799QKF.bam
mv ENCFF799QKF.bam b_cell_low.bam
# Peaks
wget https://www.encodeproject.org/files/ENCFF507JIF/@@download/ENCFF507JIF.bed.gz
gunzip ENCFF507JIF.bed.gz
mv ENCFF507JIF.bed b_cell_low.bed