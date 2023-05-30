#!/bin/bash

for i in size18-35_trimmed/*.gz
do
	h=${i##*/}
	k=${h%.fastq.gz}
	echo $k
	gunzip -c $i | bowtie -q -v 1 -p 8 -S -a -m 100 --best --strata reference/ddis_genomic - | samtools view -o mapped/${k}.bam
done