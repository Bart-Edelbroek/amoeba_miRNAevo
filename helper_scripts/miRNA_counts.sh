featureCounts -a ../data_in/ddis/sRNAs_ddis_mature.gff -T 8 -O -M --fraction --fracOverlap 1  --fracOverlapFeature 0.9 -t miRNA-5p -g Name -o ../data_in/miRNA_5p_counts.txt mapped/* 
featureCounts -a ../data_in/ddis/sRNAs_ddis_mature.gff -T 8 -O -M --fraction --fracOverlap 1  --fracOverlapFeature 0.9 -t miRNA-3p -g Name -o ../data_in/miRNA_3p_counts.txt mapped/* 