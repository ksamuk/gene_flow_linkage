fit_linear_model_dxy <- function (matched.all){
	
	if (sum(!is.na(matched.all$dxy.outlier)) < 100){
		warning("Low data, aborting model fitting")
	} else{
		model <- glm(as.numeric(dxy.outlier) ~ recomb_rate + ds + gene_count,
								 na.action = "na.omit",
								 family = binomial,
								 data = matched.all)
		return(model)
	}
	
}