cat ../../species.txt | while read i
do
cmsearch --noali --tblout align_all.txt --oskip amoebozoa.cm ../${i}/*_genomic.fna
grep '!' align_all.txt | grep '+' | awk -v OFS='\t' '{print $1,$8,$9,$3,$16,$10}' > ../${i}/ncRNA.bed
grep '!' align_all.txt | grep -v '+' | awk -v OFS='\t' '{print $1,$9,$8,$3,$16,$10}' >> ../${i}/ncRNA.bed
Rscript rfam_to_bed.R ../${i}/ncRNA.bed $i
done