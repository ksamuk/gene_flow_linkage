fit_linear_model_fst <- function (matched.all){
	
	model <- glm(as.numeric(fst.outlier) ~ recomb_rate + ds + gene_count,
							 na.action = "na.omit",
							 family = binomial,
							 data = matched.all)
	
	return(model)
}