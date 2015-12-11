########################################
# Calculate a two-sided monte carlo
# style p-value from a permuted distribution
# and an empirical mean
########################################


two_side_p <- function(dist, mean){
	p1 <- (sum(dist > mean)+1) / (length(dist)+1)
	p2 <- (sum(dist < mean)+1) / (length(dist)+1)
	p <- min(p1, p2)*2
	return(p)
}