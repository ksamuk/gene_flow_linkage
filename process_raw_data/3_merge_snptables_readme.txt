This is a description of how the accompanying scripts are used to merge multiple large tab separated snp tables into one large table that is filtered so that each site must be represented by more than one population.

#All the sites in all the tables are collected and ordered.
usage: perl ./3a_make_master_site_list.pl list_of_snp_tables.txt > all_sites.txt

#Each individual snp table is remastered to include all sites. This puts large amounts of missing data.
usage: cat all_sites.txt | perl ./3b_merge_snp_tables.v2.pl snp_table_1.tab > snp_table_1_remastered.tab

#After remastering, snp tables are merged using command line tools
usage: paste snp_table_1_remastered.tab <(cut -f 3- snp_table_2_remastered.tab ) > snp_table_1.2_remastered.tab

#Final SNPtable is filtered so each site is in more than one population.
usage: cat snp_table_all_remastered.tab | perl ./3c_filter_large_snptable.pl population_file.txt > snp_table_all_filtered.tab
