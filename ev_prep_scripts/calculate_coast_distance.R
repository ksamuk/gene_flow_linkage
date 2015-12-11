############################################################
# calculate euclidian & 'least cost' distance between populations
# KS Aug 2015
############################################################

rm(list =ls())

################################################################################
# Libraries and functions
################################################################################

library(parallel)
library(marmap)
library(dplyr)
library(fossil)

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select

################################################################################
# Input files
################################################################################

pop.dat <- read.csv("meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
stats.file.names <- list.files("stats/75k_all")

################################################################################
# Initialize marmap bathy data 
################################################################################

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

# preload and transform bathy data
# this is a big bottleneck if placed in the function
bat1 <- getNOAA_bathy_prefix(0, -45, 30, 80, res = 8, keep = TRUE, antimeridian = TRUE, prefix = "meta_data/")
bat2 <- getNOAA_bathy_prefix(180, -180, 30, 80, res = 8, keep = TRUE, prefix = "meta_data/")

trans.1.file <- "meta_data/bathy.1.trans"
trans.2.file <- "meta_data/bathy.2.trans"

# if the transformed files don't exist, create them and save them to disk
if(!file.exists(trans.1.file)){
	tr1 <- trans.mat(bat1, min.depth = -10, max.depth = -800)
	save(tr1, file = trans.1.file)
} else{
	load(trans.1.file )
}

if(!file.exists(trans.2.file)){
	tr2 <- trans.mat(bat2, min.depth = -10, max.depth = -800)
	save(tr2, file = trans.2.file)
} else{
	load(trans.2.file )
}

################################################################################
# Run core function (compute_lc_distance_stats_file)
################################################################################

# run the above function for all files, on three cores 
# might not work in *NIX?

cl <- makeCluster(getOption("cl.cores", 3))

clusterEvalQ(cl, {
	library(parallel)
	library(marmap)
	library(dplyr)
	library(fossil)
	list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible
	pop.dat <- read.csv("meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
	stats.file.names <- list.files("stats/75k_all")
	blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
	greys <- c(grey(0.6), grey(0.93), grey(0.99))
})

distances.df <- parLapply(cl = cl, stats.file.names, compute_lc_distance_stats_file, 
													bathy1 = bat1, bathy2 = bat2, trans1 = tr1, trans2 = tr2,
													plot.map = TRUE)

stopCluster(cl)

# single core version
# works for all OSs
#distances.df <- lapply(stats.file.names, compute_lc_distance_stats_file, 
#											 bathy1 = bat1, bathy2 = bat2, trans1 = tr1, trans2 = tr2,
#											 plot.map = TRUE)

# clean up results and write to file
distances.df <- bind_rows(distances.df)

distances.df <- distances.df %>%
	unique

write.table(distances.df, file = "meta_data/pop_geo_distances.txt", 
						quote = FALSE, row.names = FALSE)
