fit_linear_model_dxy_quant <- function (matched.all){
	
	model <- lm(dxy ~ recomb_rate + ds + gene_count,
							 na.action = "na.omit",
							 data = matched.all)
	
	return(model)
}