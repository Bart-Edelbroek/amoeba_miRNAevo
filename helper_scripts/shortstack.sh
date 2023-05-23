#!/bin/bash


cat species.txt | while read i
do
	echo "processing" $i
	Shortstack --readfile size18-35_trimmed/$i*.fastq.gz --outdir shstack_out/$i/ --pad 100 \
	--bowtie_cores 8 --mismatches 1 -mmap f --genomefile references/$i/*_comb.fa
 
done