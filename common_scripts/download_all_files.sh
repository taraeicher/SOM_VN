#GM12878 - https://www.encodeproject.org/experiments/ENCSR000EJD/
#BAM 
wget https://www.encodeproject.org/files/ENCFF577DUO/@@download/ENCFF577DUO.bam
cp ENCFF577DUO.bam GM12878.bam
#Peak
wget https://www.encodeproject.org/files/ENCFF588OCA/@@download/ENCFF588OCA.bed.gz
gunzip ENCFF588OCA.bed.gz
cp ENCFF588OCA.bed GM12878_peaks.bed
#ChromHMM mnemonics
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E116_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E116_15_coreMarks_hg38lift_mnemonics.bed.gz
cp E116_15_coreMarks_hg38lift_mnemonics.bed GM12878_mnemonics.bed

#A549 - https://www.encodeproject.org/experiments/ENCSR136DNA/
#BAM
wget https://www.encodeproject.org/files/ENCFF473YHH/@@download/ENCFF473YHH.bam
wget https://www.encodeproject.org/files/ENCFF809KIH/@@download/ENCFF809KIH.bam
wget https://www.encodeproject.org/files/ENCFF961WXW/@@download/ENCFF961WXW.bam
wget https://www.encodeproject.org/files/ENCFF821UUL/@@download/ENCFF821UUL.bam
bamtools merge -in ENCFF473YHH.bam -in ENCFF809KIH.bam -in ENCFF961WXW.bam -in ENCFF821UUL.bam -out A549.bam
#Peak
wget https://www.encodeproject.org/files/ENCFF135JRM/@@download/ENCFF135JRM.bed.gz
gunzip ENCFF135JRM.bed.gz
wget https://www.encodeproject.org/files/ENCFF698UAH/@@download/ENCFF698UAH.bed.gz
gunzip ENCFF698UAH.bed.gz
wget https://www.encodeproject.org/files/ENCFF079DJV/@@download/ENCFF079DJV.bed.gz
gunzip ENCFF079DJV.bed.gz
wget https://www.encodeproject.org/files/ENCFF045PYX/@@download/ENCFF045PYX.bed.gz
gunzip ENCFF045PYX.bed.gz
bedtools intersect -a ENCFF135JRM.bed -b ENCFF698UAH.bed ENCFF079DJV.bed ENCFF045PYX.bed > A549_peaks.bed
#ChromHMM mnemonics
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E114_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E114_15_coreMarks_hg38lift_mnemonics.bed.gz
cp E114_15_coreMarks_hg38lift_mnemonics.bed A549_mnemonics.bed

#H1-hESC - https://www.encodeproject.org/experiments/ENCFF546PJU/
#BAM
wget https://www.encodeproject.org/files/ENCFF546PJU/@@download/ENCFF546PJU.bam
cp ENCFF546PJU.bam H1.bam
#Peak
wget https://www.encodeproject.org/files/ENCFF030XPN/@@download/ENCFF030XPN.bed.gz
gunzip ENCFF030XPN.bed.gz
cp ENCFF030XPN.bed H1_peaks.bed
#ChromHMM mnemonics
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E003_15_coreMarks_hg38lift_mnemonics.bed.gz
cp E003_15_coreMarks_hg38lift_mnemonics.bed. H1_mnemonics.bed

#Brain tissue - https://www.encodeproject.org/experiments/ENCSR649KBB/
#BAM
wget https://www.encodeproject.org/files/ENCFF002QDM/@@download/ENCFF002QDM.bam
cp ENCFF002QDM.bam Brain.bam
#Peak
wget https://www.encodeproject.org/files/ENCFF286XIL/@@download/ENCFF286XIL.bed.gz
gunzip ENCFF286XIL.bed.gz
cp ENCFF286XIL.bed Brain_peaks.bed
#ChromHMM mnemonics
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E081_15_coreMarks_hg38lift_mnemonics.bed.gz
gunzip E081_15_coreMarks_hg38lift_mnemonics.bed.gz
cp E081_15_coreMarks_hg38lift_mnemonics.bed Brain_mnemonics.bed