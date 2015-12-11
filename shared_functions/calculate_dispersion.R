calculate_dispersion <- function (lg) {
	
	dispersion <- lg %>%
		select(gen.pos) %>% 
		arrange(gen.pos) %>%
		unlist %>% 
		diff %>% 
		(function(x) return(var(x,na.rm = TRUE)/mean(x,na.rm = TRUE)))
	
	return(dispersion)
	
}