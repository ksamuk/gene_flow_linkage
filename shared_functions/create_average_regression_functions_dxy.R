create_average_regression_functions_dxy <- function(coeff.dat){
	make_regression_function <- function(intercept, slope){
		
		reg_function <- function(x){
			
			ylog <- intercept + slope*x 
			return(1 - 1/(1 + exp(ylog)))
			
		}
		return(reg_function)
	}
	
	avg_regression_coeffs <- coeff.dat %>%
		filter(!is.na(recomb_rate_dxy)) %>%
		group_by(group2.new) %>%
		summarise(mean_recomb_rate = mean(recomb_rate_dxy),
							mean_intercept = mean(intercept_dxy))
	
	avg_regression_functions <- list()
	
	avg_regression_functions[[1]] <- make_regression_function(avg_regression_coeffs$mean_intercept[1], avg_regression_coeffs$mean_recomb_rate[1])
	avg_regression_functions[[2]] <- make_regression_function(avg_regression_coeffs$mean_intercept[2], avg_regression_coeffs$mean_recomb_rate[2])
	avg_regression_functions[[3]] <- make_regression_function(avg_regression_coeffs$mean_intercept[3], avg_regression_coeffs$mean_recomb_rate[3])
	avg_regression_functions[[4]] <- make_regression_function(avg_regression_coeffs$mean_intercept[4], avg_regression_coeffs$mean_recomb_rate[4])
	
	names(avg_regression_functions) <- avg_regression_coeffs$group2.new
	return(avg_regression_functions)
}

