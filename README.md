# amoeba_miRNAevo
Analysis of miRNAs in amoebozoa, and their characteristics

## microRNA curation
Potential miRNA clusters, generated using Shortstack (https://github.com/MikeAxtell/ShortStack), are analyzed with **miRNA_curation.py**. An example run of the miRNA curation can be performed by running the 'example_run/explore_all.sh' and 'example_run/refine_all.sh' scripts in order. A sample of the sRNA sequences from *D. discoideum* is included in the 'example/example_sizezselected_data' folder. A sample of miRNA candidates are included in 'example_run/pot-miRNAs_merged.fa' to facilitate the example run. <br />
Prerequisites: Python 3, Biopython, RNAlib, matplotlib, numpy

## microRNA analysis
Following miRNA curation, high confidence miRNAs are identified and further analysed in **amoeba_miRNAevo.Rmd**. The data for the analysis is located in the 'data_in' and is mostly derived from the miRNA_curation.py script and other scripts located in the 'helper_scripts' folder. The plots generated in the analysis is located in the 'plots' folder. Tables and other data generated from the analysis are located in the 'data_out' folder. 
