filter_data_call_outliers_function <- function (matched.all){
	
	# apply filters to data
	# these filter out extreme/bioligically implausible values
	# also throws out the sex chromosome (lg 19)
	matched.all <- matched.all %>%
		filter(!is.na(fst)) %>%
		filter(!is.infinite(fst)) %>%
		filter(recomb_rate <= 25) %>% 
		filter(var.sites >= 2) %>%
		filter(ds <= 3) %>%
		filter(gene_count <= 15) %>%
		filter(lg!=19) %>%
		mutate(dxy_adj = ifelse(sites >= 500, dxy, NA))
	
	# force negative fst to be zero
	matched.all$fst[matched.all$fst < 0] <- 0 
	
	# compute outliers
	matched.all <- matched.all %>%
		mutate(fst.outlier = is.outlier(fst))%>%
		mutate(dxy.outlier = is.outlier(dxy_adj)) %>%
		mutate(both.outlier = fst.outlier & dxy.outlier)%>%
		mutate(hs = (hexp1+hexp2)/2)
	
	return(matched.all)
}