save_pvals <- function(permutation_output, stat, group_type){
	
	pvals <- permutation_output$pvals
	pvals$stat <- stat
	pvals$group_type <- group_type
	pvals$observed.means <- unlist(permutation_output$observed.means[,2])
	pvals$permuted.means <- unlist(permutation_output$permuted.means[,2])
	
	return(pvals)
	
}
