wget -O TE_all.fa "https://www.girinst.org/protected/repbase_extract.php?division=Eukaryota&customdivision=&rank=&type=Transposable+Element&autonomous=1&nonautonomous=1&simple=1&format=FASTA&sa=Download"
grep -E "Dictyostelium discoideum|Dictyostelium firmibasis|Dictyostelium lacteum|Dictyostelium fasciculatum|Polysphondylium pallidum|Acanthamoeba lenticulata|Acanthamoeba castellanii|Acytostelium subglobosum|Physarum polycephalum" TE_all.fa | cut -c 2- > amoebozoa_IDs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' amoebozoa_IDs.txt TE_all.fa > TE_amoebozoa.fa
awk 'BEGIN{RS=">";ORS=""}length($0)>500{print ">"$0}' TE_amoebozoa.fa > TE_amoebozoa_500.fa
rm TE_all.fa