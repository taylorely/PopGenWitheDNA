# PopGenWitheDNA
Scripts and data for my first chapter comparing eDNA methods to traditional tissue-based methods for population genetic and genetic diversity analyses for marine vertebrate communities in the Main Hawaiian Islands.
This paper can be found:
add later

The file processing_eDNA_data.R has the code for loading and reformating JAMP, ecotag, and blast outputs. It also contains all the code to create the proxies for individuals and run sequence diversity statistics, neutrality statistics, and population structure statistics. Since I was repeating the same functions for multiple datasets (different filter quality thresholds and different proxies), I created functions in order to easily repeat these analyses for the multiple datasets.

The file processing_tissue_data.R uses the same functions as for the eDNA file but with the data from previous studies. This data was downloaded from genbank, dryad, or other online repositories. All data was reformatted to a fasta file. An example file for the full original length and the shortened length can be found in the input_files folder (Abuabd_long.fasta and Abuabd_short.fasta). The format for the names in the fasta file is population_haplotypeid_haplotypecount.

The file HaplotypeMapsPieCharts.R has the code for making haplotype networks, haplotype pie charts, and stacked barplots of haplotype frequencies using the same colors for the same haplotype detected with both eDNA and tissue-based methods. These figures were used for Fig. 1, 2, S1, and S2

The file Correlations_PairedTtests.R has the code to run Spearman rank correlations and paired t-tests used for Table 2. It also makes the figures for Fig. 3 and 4 and S3-S13.

The folder files_for_comparisons contains the input files for Correlations_PairedTtests.R
The folder input_files has the input files for processing eDNA_data.R and example input files for processing_tissue_data.R
The folder example outputs are example outputs from processing_eDNA_data.R