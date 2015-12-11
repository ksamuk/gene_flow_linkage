nuke_duplicate_cluster_files <- function() {
	
	out.folder <- "analysis_ready/clustering_fst_new"
	out.files <- list.files(out.folder)
	
	# check which SNP files (stats.folder) still need to be processed
	out.files.full <- list.files(out.folder, full.names = TRUE)
	out.files <- list.files(out.folder)
	out.files.exist <- out.files %>% gsub("\\d","",.) %>% gsub("-","",.) %>% gsub(".gz_.clustered.txt",".txt.gz",.)
	
	# nuke duplicates
	out_df <- data.frame(full_name = out.files.full)
	out_df$short_name <- out.files.exist
	out_df$date_slug <- lapply(strsplit(out.files, ".gz_"), function(x)x[2])
	
	dup_log <- duplicated(out_df$short_name)
	print(paste(sum(dup_log), "duplicates found. removing..."))
	
	dups <- as.character(out_df$full_name[duplicated(out_df$short_name)])
	
	lapply(dups, file.remove) %>% invisible

}
