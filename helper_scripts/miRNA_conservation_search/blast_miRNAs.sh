for i in *.fa
do
	s=${i##*/}
	p=${s%.*}
	makeblastdb -in ${p}.fa -dbtype nucl -out blast/$p
	blastn -task blastn-short -db blast/$p -query ../miRNAs_withstar.fa -outfmt "6 sseqid sstart send qseqid length evalue bitscore" -evalue 2 -out txt/${p}.txt
	Rscript ../output_to_bed.R txt/${p}.txt bed/${p}.bed
done