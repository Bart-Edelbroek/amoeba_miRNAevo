cat species.txt | while read i
do
	echo $i >> "length_profiles/mapping_percentages.txt"
	samtools view shstack_out_fraction/${i}/merged_alignments.bam | grep -E 'contig_142|contig_149|contig_158' | wc -l >> "length_profiles/mapping_percentages.txt"
	samtools view shstack_out_fraction/${i}/merged_alignments.bam | grep -E 'XY:Z:N|XY:Z:M|XY:Z:O' | wc -l >> "length_profiles/mapping_percentages.txt"
	samtools view shstack_out_fraction/${i}/merged_alignments.bam | grep -vE 'contig_142|contig_149|contig_158|XY:Z:N|XY:Z:M|XY:Z:O' | wc -l >> "length_profiles/mapping_percentages.txt"
done