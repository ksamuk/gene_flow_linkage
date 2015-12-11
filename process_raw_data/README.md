# Processing raw reads
These scripts are designed to start with raw fastq reads and ends with Fst and Dxy scores for each site and genomic window for all population comparisons. 

###1a_process_fastq.sh
-This script takes the raw reads, aligns them to the stickleback genome using BWA and stampy, processes the BAM files using samtools and GATK, and calls SNPs using GATK. 

###1b_GBS_fastq_demultiplexer.pl
-This script is part of the 1a_process_fastq.sh pipeline and splits raw GBS reads into individuals as well as splits off barcodes and adapter sequences.

###2_vcf2vertical_dep_GATK33_full.pl
-This script takes a vcf file and converts it to a flat tab separated format, including invariant sites. It requires 5 reads per genotype and removes indels and tri/quad-allelic sites.

###3_merge_snptables_readme.txt
-After a .tab file is created for each project, the SNP tables need to be merged. This is complicated by the differences in coverage across the genome. This readme file explains the process of merging the snp tables.

###3a_make_master_site_list.pl
-This script creates a master list of all possible sites with data

###3b_merge_snp_tables.v2.pl
-This script takes a list of possible sites and remasters a snp table so that all those sites are included (i.e. creates rows of missing data so each snp table has equal rows)

###3c_filter_large_snptable.pl
-This script filters a large snp table to make sure there are genotypes in each site in multiple populations (i.e. to remove data that won't be used)

###4_Run_all_meta_comparisons.pl
-This script prepares scripts to run all possible pairwise comparisons to calculate Fst/dxy in windows. It is designed to be run on WestGrid servers.

###4a_SNPtable2Fstdxy_filtered_speedy.pl
-This script calculates weir and cockerhams Fst, as well as Dxy on each site.

###4b_slidingwindow_v5.pl
-This script runs a sliding window across the genome and averages Fst and Dxy in each window.

###grouplist_example.txt
-This file lists the populations along with location and ecotype

###samplelist_example.txt
-This file lists the sample names and population.

###snptable_example.tab
-This is an example of the tab separated snp table.
