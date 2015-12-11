add_region_data <- function(coeff.dat){

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
	
	#make new (relaxed) groups 
	coeff.dat$group2 <- paste0(coeff.dat$geography2,"_",coeff.dat$ecology)
	
	return(coeff.dat)

}
