# make supplemental datafiles 1, 2, 3

# read in raw data
dat.fst <- read.delim(file = "analysis_ready/75k_stats_model_fits_fst.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
dat.dxy <- read.delim(file = "analysis_ready/75k_stats_model_fits_dxy.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
dat.hs <- read.delim(file = "analysis_ready/75k_stats_model_fits_hs.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")


# format for nice output
dat.fst <- with(dat.fst, data.frame(dat.fst[,1:4],
																		region1 = reg1,
																		region2 = reg2,
												            geography = geography2, 
																		group = group2,
																		dat.fst[,7:10]))

dat.dxy <- with(dat.dxy, data.frame(dat.dxy[,1:4],
																		region1 = reg1,
																		region2 = reg2,
																		geography = geography2, 
																		group = group2,
																		dat.dxy[,7:10]))

dat.hs <- with(dat.hs, data.frame(dat.hs[,1:4],
																			 region1 = reg1,
																			 region2 = reg2,
																			 geography = geography2, 
																			 group = group2,
																			 dat.hs[,7:10]))

# write to file
write.table(dat.fst, file = "figures/supplemental_data_1_model_fits_fst.txt", quote = FALSE, row.names = FALSE)
write.table(dat.dxy, file = "figures/supplemental_data_2_model_fits_dxy.txt", quote = FALSE, row.names = FALSE)
write.table(dat.hs, file = "figures/supplemental_data_3_model_fits_hs.txt", quote = FALSE, row.names = FALSE)

