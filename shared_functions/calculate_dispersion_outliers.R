calculate_dispersion_outliers <- function (lg) {
	
	start <- min(lg$gen.pos, na.rm = TRUE)
	end <- max(lg$gen.pos, na.rm = TRUE)
	
	
	outlier.positions <- lg %>% 
		filter(fst.outlier == TRUE) %>%
		select(gen.pos) %>%
		unlist
	outlier.positions <- c(start, outlier.positions, end)
	
	outlier.dispersion <- outlier.positions %>%
		sort %>%
		diff %>% 
		(function(x) return(var(x,na.rm = TRUE)/mean(x,na.rm = TRUE)))
	
	return(outlier.dispersion)
	
}