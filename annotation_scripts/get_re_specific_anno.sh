# Load modules.
	module load cuDNN/7.6.5/CUDA-10.1 python/3.7
	module load bedtools
	

#for sample in "b_cell_low" "b_cell_high" "Brain_low" "H1_high" "HeLa_low" "HeLa_high" "heart_low" "heart_high" "stomach_low" "stomach_high";
for sample in "GM12878" "Brain_high" "H1_low" "A549";
do

#Variables
	for i in {1..22};
	do
		awk '($4 == "Enhancer")' /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno$i.bed > /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno${i}_known_enhancer.bed
		awk '($4 == "Promoter")' /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno$i.bed > /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno${i}_known_promoter.bed
		awk '($4 == "Weak")' /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno$i.bed > /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno${i}_known_weak.bed
	done
	cat /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno*_known_enhancer.bed > /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno_all_known_enhancer.bed
	cat /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno*_known_promoter.bed > /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno_all_known_promoter.bed
	cat /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno*_known_weak.bed > /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno_all_known_weak.bed
	bedtools intersect -a /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno_all_known_enhancer.bed -b /data/eichertd/som_vn_data/peaks/$sample.bed -u > annotations_enhancer_${sample}_peak.bed
	bedtools intersect -a /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno_all_known_promoter.bed -b /data/eichertd/som_vn_data/peaks/$sample.bed -u > annotations_promoter_${sample}_peak.bed
	bedtools intersect -a /data/eichertd/som_vn_data/annotated_consolidated_$sample/anno_all_known_weak.bed -b /data/eichertd/som_vn_data/peaks/$sample.bed -u > annotations_weak_${sample}_peak.bed
done 