cat ../../species.txt | while read i
do
	makeblastdb -in ../$i/${i}*_genomic.fna -dbtype nucl
	tblastx -query TE_amoebozoa_500.fa -db ../$i/${i}*_genomic.fna -outfmt "6 sseqid sstart send qseqid length evalue bitscore" -evalue 10E-15 > ${i}_out.txt
	Rscript output_to_bed.R ${i}_out.txt ../$i/annot_TEmerged.bed
done