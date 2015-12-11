fit_linear_model_fst_quant <- function (matched.all){
	
	model <- lm(fst ~ recomb_rate + ds + gene_count,
							 na.action = "na.omit",
							 data = matched.all)
	
	return(model)
}