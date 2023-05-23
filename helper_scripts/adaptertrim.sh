#!/bin/bash

for i in merged_fastq/*.gz
do
	cutadapt -a TGGAATTCTCGGGTGCCAAGG -j 0 -o size18-35_trimmed/${i##*/} -m 18 -M 35 $i
	echo size18-35_trimmed/${i##*/}
done