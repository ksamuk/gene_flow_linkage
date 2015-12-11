initialize_coeff_dat_files<- function(){
	
	rename <- dplyr::rename

	# choose either fst or fst/dxy 
	coeff.dat.fst <- read.table(file = "analysis_ready/75k_stats_model_fits_fst.txt", header = TRUE, stringsAsFactors = FALSE, sep="\t")
	coeff.dat.dxy <- read.table(file = "analysis_ready/75k_stats_model_fits_dxy.txt", header = TRUE, stringsAsFactors = FALSE, sep="\t")
	coeff.dat.hs <- read.table(file = "analysis_ready/75k_stats_model_fits_hs.txt", header = TRUE, stringsAsFactors = FALSE, sep="\t")
	
	coeff.dat.fst <- coeff.dat.fst %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		rename(recomb_rate_fst = recomb_rate) %>%
		rename(intercept_fst = intercept) %>%
		rename(model_warning_fst = linear_model_warning) %>%
		select(comparison, pop1, ecotype1, pop2, ecotype2, 
					 reg1, reg2, geography, geography2, ecology, 
					 group, group2, intercept_fst, recomb_rate_fst, n_windows_fst, model_warning_fst)
	
	coeff.dat.dxy <- coeff.dat.dxy %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		rename(recomb_rate_dxy = recomb_rate) %>%
		rename(intercept_dxy = intercept) %>%
		rename(model_warning_dxy = linear_model_warning) %>%
		select(comparison, intercept_dxy, recomb_rate_dxy, n_windows_dxy, model_warning_dxy)
	
	coeff.dat.hs <- coeff.dat.hs %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		rename(recomb_rate_hs = recomb_rate) %>%
		rename(intercept_hs = intercept) %>%
		rename(model_warning_hs = linear_model_warning) %>%
		select(comparison, intercept_hs, recomb_rate_hs, model_warning_hs)
	
	
	coeff.dat <- left_join(coeff.dat.fst, coeff.dat.dxy)
	coeff.dat <- left_join(coeff.dat, coeff.dat.hs)
	
	## add in region data (for looser geography)
	region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)
	
	# make short pop codes
	region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist
	
	# associate pops in coeff dat with regions in region dat
	region.sub <- region.dat[,c(1,4)]
	names(region.sub) <- c("pop1","reg1")
	coeff.dat$reg1 <- region.sub$reg1[match(coeff.dat$pop1, region.sub$pop1)]
	
	region.sub <- region.dat[,c(1,4)]
	names(region.sub) <- c("pop2","reg2")
	coeff.dat$reg2 <- region.sub$reg2[match(coeff.dat$pop2, region.sub$pop2)]
	
	# make new geographic categories
	coeff.dat$geography2 <- ifelse(coeff.dat$reg1==coeff.dat$reg2, "para", "allo")
	
	#make new groups :o
	coeff.dat$group2 <- paste0(coeff.dat$geography2,"_",coeff.dat$ecology)
	
	group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
	group.rename <- c("Allopatry\nDivergent", "Allopatry\nParallel", "Gene Flow\nDivergent", "Gene Flow\nParallel")
	coeff.dat$group2.new <- group.rename[match(coeff.dat$group2, group.old.names)]
	coeff.dat$group.new <- group.rename[match(coeff.dat$group, group.old.names)]
	
	# filter out comparisons that had fewer than 200 windows with > 500 sites (n_windows_dxy is a count of these)
	
	return(coeff.dat)
}