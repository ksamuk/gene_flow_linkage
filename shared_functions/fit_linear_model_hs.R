fit_linear_model_hs <- function (matched.all){
	
		model <- lm(hs ~ recomb_rate + ds + gene_count,
								 na.action = "na.omit",
								 data = matched.all)
		return(model)
	
}