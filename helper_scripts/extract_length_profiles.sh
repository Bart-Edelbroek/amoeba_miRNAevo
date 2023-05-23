extract_lengths(){
	> references/${1}/annot_mRNA.bed
	> references/${1}/annot_intergenicRNA.bed
	cat references/${1}/annot_*.bed > references/${1}/merged.bed
	bedtools sort -g references/${1}/chrom.list -i references/${1}/merged.bed > references/${1}/merged_sorted.bed
	bedtools subtract -a references/${1}/mRNA.bed -b references/${1}/merged_sorted.bed | bed12toBed6 > references/${1}/annot_mRNA.bed
	cat references/${1}/merged_sorted.bed references/${1}/annot_mRNA.bed | bedtools sort -g references/${1}/chrom.list -i - | bedtools complement -i - -g references/${1}/${1}_bed.genome | grep -vE 'contig_142|contig_149|contig_158' > references/${1}/annot_intergenicRNA.bed
	> length_profiles/$1/length_profile.txt
	c="complementary_"
	for m in references/${1}/annot_*.bed
	do
		l=${m##*_}
		s=${l%.*}
		echo "${1}_$s"
		if [[ "$s" == "intergenicRNA" ]]; then
			bedtools intersect -f 0.5 -b $m -a shstack_out/${1}/merged_alignments.bam | samtools view | grep -v '^@' | awk '{print length($10), 1}' | awk -v b="$s" 'NR == 0 { print; next } { a[$1] += $2 } END { for (i in a) { print  a[i], i, b; } } ' | sort -k 2 >> length_profiles/$1/length_profile.txt
			continue
		fi
		bedtools intersect -f 0.5 -s -b $m -a shstack_out/${1}/merged_alignments.bam | samtools view | grep -v '^@' | awk '{print length($10), 1}' | awk -v b="$s" 'NR == 0 { print; next } { a[$1] += $2 } END { for (i in a) { print  a[i], i, b; } } ' | sort -k 2 >> length_profiles/$1/length_profile.txt
		q="$c$s"
		echo "${1}_$q"
		bedtools intersect -f 0.5 -S -b $m -a shstack_out/${1}/merged_alignments.bam | samtools view | grep -v '^@' | awk '{print length($10), 1}' | awk -v b="$q" 'NR == 0 { print; next } { a[$1] += $2 } END { for (i in a) { print  a[i], i, b; } } ' | sort -k 2 >> length_profiles/$1/length_profile.txt
	done
}
export -f extract_lengths

parallel -u -a species.txt extract_lengths