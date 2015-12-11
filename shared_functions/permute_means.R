################################################################################
# shuffle a statistic over 'groups' defined in the data, and return the mean
# statistic for each group (i.e. permutation)
################################################################################

permute_means <- function(data, stat, group_type) {
	# shuffle groups
	data[,group_type] <- sample(data[,group_type] %>% unlist, length(data[,group_type] %>% unlist))
	
	# calculate shuffled group means
	mean <- data %>% 
		group_by_(interp(group_type)) %>% 
		summarise_(mean_stat = interp(~ mean(var, na.rm = TRUE), var = as.name(stat))) %>% 
		ungroup
	names(mean)[2] <- stat
	return(mean)
}


