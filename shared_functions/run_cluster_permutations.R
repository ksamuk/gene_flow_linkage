run_cluster_permutations <- function(cluster.df, stat, group_type, n_permutations, collapsed = FALSE){
	
	if(collapsed == FALSE){
	#crunch multi-lg metrics into single measurements
	cluster.df <- cluster.df %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		group_by_(interp(group_type), ~comparison) %>%
		summarise_(mean_stat = interp(~ mean(var, na.rm = TRUE), var = as.name(stat))) %>%
		ungroup
		names(cluster.df)[3] <- stat
	}
	
	# run the permutation function above for 10,000 interations and bind into df
	permuted.means.list <- replicate(n_permutations, 
																	 permute_means(cluster.df, stat, group_type), simplify = FALSE)
	permuted.means.df <- bind_rows(permuted.means.list)
	
	#calculate observed means
	observed.means <- cluster.df %>% 
		group_by_(interp(group_type)) %>%
		summarise_(mean_stat = interp(~ mean(var, na.rm = TRUE), var = as.name(stat))) %>% 
		ungroup
	
	permuted.means <- permuted.means.df %>% 
		group_by_(interp(group_type)) %>%
		summarise_(mean_stat = interp(~ mean(var, na.rm = TRUE), var = as.name(stat))) %>% 
		ungroup
	
	# ecdfs for plotting
	ecdf.df <- permuted.means.df %>% 
		group_by_(interp(group_type)) %>% 
		do_(ecdf = interp(~ ecdf(.$var), var = as.name(stat)))
	
	ecdf.df$obs <- observed.means[,"mean_stat"] %>% unlist
	
	# calculate pvalues (used in paper)
	pvals <- list()
	group_names <- unlist(unique(permuted.means.df[,1]))
	for (i in 1:length(group_names)){
		df <- subset(permuted.means.df, permuted.means.df[,1] == group_names[i])
		pvals[[i]] <- data.frame(pvalue = two_side_p(unlist(df[,2]), observed.means$mean_stat[i]), group = group_names[i])
	}

	output <- list(permuted.means.df = permuted.means.df, 
								 observed.means = observed.means,
								 permuted.means = permuted.means,
								 ecdf.df = ecdf.df, 
								 pvals = bind_rows(pvals))
	return(output)
}
