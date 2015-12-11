make_regression_function <- function(intercept, slope){
	
	reg_function <- function(x){
		
		ylog <- intercept + slope*x 
		return(1 - 1/(1 + exp(ylog)))
		
	}
	return(reg_function)
}