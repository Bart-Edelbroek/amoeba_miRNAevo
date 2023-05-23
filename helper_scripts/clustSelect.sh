#!/bin/bash
folder=shstack_out
cat species.txt | while read i
do
	Rscript clustSelect.R $folder/$i/Results.txt $folder/$i/ShortStack_All.gff3 25 $folder/$i/pot_clusters_overhang.gff3
	Rscript clustSelect.R $folder/$i/Results.txt $folder/$i/ShortStack_All.gff3 0 $folder/$i/pot_clusters.gff3
	gff2bed < $folder/$i/pot_clusters_overhang.gff3 > $folder/$i/pot_clusters_overhang.bed
	gff2bed < $folder/$i/pot_clusters.gff3 > $folder/$i/pot_clusters.bed
	echo $i
	bedtools getfasta -fi references/$i/*_comb.fa -fo $folder/$i/pot-miRNAs_overhang.fa -bed $folder/$i/pot_clusters_overhang.bed -nameOnly -s 
	bedtools getfasta -fi references/$i/*_comb.fa -fo $folder/$i/pot-miRNAs.fa -bed $folder/$i/pot_clusters.bed -nameOnly -s 
	mkdir -p miranalysis/${i}
	awk '/^>/{p=seen[$0]++}!p' $folder/$i/pot-miRNAs_overhang.fa $folder/$i/pot-miRNAs.fa > miranalysis/$i/pot-miRNAs_merged.fa
done
