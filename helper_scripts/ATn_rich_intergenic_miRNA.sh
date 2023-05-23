echo "species,type,A,T,C,G,N,total" > ATn_rich_intergenic_miRNA.csv
cat ../species.txt | while read i
do
	bedtools getfasta -fi $i/$i*_genomic.fna -bed $i/annot_intergenicRNA.bed | grep -v '^>' > nucleotides_intergenic.txt
	bedtools getfasta -fi $i/$i*_genomic.fna -bed $i/annot_mRNA.bed | grep -v '^>' > nucleotides_ORF.txt
	bedtools getfasta -fi $i/$i*_genomic.fna -bed $i/primir.bed | grep -v '^>' > nucleotides_miRNA.txt
	for s in nucleotides_intergenic.txt nucleotides_ORF.txt nucleotides_miRNA.txt
	do
		A=$(grep -oE 'A|a' $s | wc -l)
		T=$(grep -oE 'T|t' $s | wc -l)
		C=$(grep -oE 'C|c' $s | wc -l)
		G=$(grep -oE 'G|g' $s | wc -l)
		N=$(grep -oE 'N|n' $s | wc -l)
		total=$(expr $A + $T + $C + $G + $N)
		t=${s##*_}
		echo $i ${t%.*} $A $T $C $G $N $total | tr ' ' , >> ATn_rich_intergenic_miRNA.csv
	done
done