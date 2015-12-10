####Match EVs to 75k stats files, fit models, write model results to file

rm(list=ls())

################################################################################
# Check and parse command arguments (if any)
################################################################################

if (length(commandArgs(TRUE)) > 0){
	
	args <- commandArgs(TRUE)
	args <- as.list(args)
	
	if(!(args[[1]] %in% c("fst","dxy","fst_dxy", "fst_quant", "dxy_quant", "hs"))){
		stop("Invalid first argument (outlier type)! Valid values are \'fst\',\'dxy\' and \'fst_dxy\' e.g.: \n\'Rscript 1_process_75k_fst.R fst\'")
	}
	
	if(length(commandArgs(TRUE)) == 2){
		
		if(is.na(as.numeric(args[[2]]))){
			stop("Invalid second argument (number of cores)! Must be an integer. e.g.:\n\'Rscript 1_process_75k_fst.R fst 4\'")
		} 
		
	}else{
		stop("Too many arguments! Run like: \n\'Rscript 1_process_75k_fst.R fst (single core mode)\'\n\'Rscript 1_process_75k_fst.R fst 4 (multicore mode)\'")
	}

}

################################################################################
# Libraries and initalizing variables
################################################################################

library("IRanges")
library("dplyr")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Input file locations
################################################################################

#ev dir location and file list
ev.dir <- file.path(getwd(),"evs")
ev.files <- list.files(ev.dir, pattern = "txt",full.names = TRUE)

#stats files location and list
stats.dir <- file.path(getwd(),"stats/75k_all")
stats.files <- list.files(stats.dir,"*sliding*", full.names = TRUE)

################################################################################
# Run main function (filter, call outliers, fit model) -- see match_evs.R
################################################################################

# the default output file name
file.name <- "analysis_ready/75k_stats_model_fst_fits.txt"

# pass arguments (if any) to mclapply/lapply
if (length(args) > 0){
	
	# choose outlier type for linear model (fst or dxy)
	if (args[[1]] == "fst"){
		linear_model_function <- fit_linear_model_fst
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	} else if (args[[1]] == "dxy"){
		linear_model_function <- fit_linear_model_dxy
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	}else if (args[[1]] == "fst_dxy"){
		linear_model_function <- fit_linear_model_fst_dxy
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	}else if (args[[1]] == "hs"){
		linear_model_function <- fit_linear_model_hs
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	} else if (args[[1]] == "fst_quant"){
		linear_model_function <- fit_linear_model_fst_quant
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	} else if (args[[1]] == "dxy_quant"){
		linear_model_function <- fit_linear_model_dxy_quant
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	} else if (args[[1]] == "hs"){
		linear_model_function <- fit_linear_model_hs
		file.name <- paste0("75k_stats_model_fits_",args[[1]],".txt")
	}
	
	# check if multicore mode and run
	if(length(args) == 2){
		cores <- args[[2]]
		
		coeff.df <- mclapply(stats.files, match_evs, 
												 linear_model_function = linear_model_function, 
												 mc.cores = cores, mc.silent = FALSE, mc.preschedule = FALSE)
		} else{
			coeff.df <- lapply(stats.files, match_evs, linear_model_function = linear_model_function)
		}
}else{
	print("No arguments provided, aborting.")
 }
#bind into a data frame
coeff.dat <- do.call("rbind", coeff.df)

#remove nas
coeff.dat <- coeff.dat %>%
  mutate(group = paste0(geography,"_",ecology))

################################################################################
# Add in region data (for relaxed geograpic classification)
################################################################################

coeff.dat <- add_region_data(coeff.dat)

################################################################################
# Write to file
################################################################################

# write to file
dir.create("analysis_ready") %>% invisible
write.table(coeff.dat, file = file.path("analysis_ready",file.name), row.names = FALSE, quote = FALSE, sep = "\t")


