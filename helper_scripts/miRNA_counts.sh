featureCounts -a ../data_in/ddis/sRNAs_ddis_mature.gff -T 8 -O -M --fraction --fracOverlap 1  --fracOverlapFeature 1 -t miRNA -g Name -o ../data_in/miRNA_counts.txt mapped/* 