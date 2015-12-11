calculate_nndist_all_lg <- function (stats.file, num_permutations, trace = FALSE) {
	
	## first, build null distributions of nndists for each linkage group:
	
	nnd.stats <- list()
	
	for (j in unique(stats.file$lg)){
		if(trace){cat(paste0("LG ", j, "..."))} 
		# subset for lg j
		stats.file.lg <- stats.file %>%
			filter(stats.file$lg == j)
		
		# the number of outliers on that lg
		num.outliers <- stats.file.lg %>%
			filter(!is.na(gen.pos)) %>%
			select(fst.outlier) %>%
			unlist %>%
			sum(na.rm = TRUE)
		
		if(trace){cat(paste0(num.outliers, " outliers."))} 
		
		if (num.outliers > 1){
			
			# draw 10000 samples of num.outliers random loci, take the mean, and return the ecdf and mean
			null.mean.nnds <- replicate(num_permutations, calculate.null.nnd(stats.file.lg, num.outliers))
			
			
			# calculate the estimate mean null nndist
			null.mean <- mean(null.mean.nnds, na.rm = TRUE)
			null.ecdf <- ecdf(null.mean.nnds)
			
			# calculate the empirical nndist for real outliers
			
			site.sample <- stats.file.lg %>%
				filter(!is.na(gen.pos)) %>%
				filter(fst.outlier == TRUE) %>%
				select(gen.pos) %>%
				arrange(gen.pos) %>%
				mutate(dist.1 = c(NA,diff(gen.pos))) %>%
				mutate(dist.2 = c(diff(sort(gen.pos)),NA))
			
			nn.dist <- rep(NA, length(site.sample$genpos))
			for (k in 1:length(site.sample$gen.pos)){
				
				if(!is.na(site.sample$dist.1[k]) & !is.na(site.sample$dist.2[k])){
					nn.dist[k] <- min(c(site.sample$dist.1[k],site.sample$dist.2[k]))
				}else if(is.na(site.sample$dist.1[k])){
					nn.dist[k] <- site.sample$dist.2[k]
				} else if(is.na(site.sample$dist.2[k])){
					nn.dist[k] <- site.sample$dist.1[k]
				}
			}
			empirical.mean.nnd <- mean(nn.dist, na.rm = TRUE)
			
			#number of total loci
			n.sites <- stats.file.lg %>% filter(!is.na(gen.pos)) %>% select(gen.pos) %>% unlist %>% length
			
			nnd.stats[[j]] <- data.frame(lg = unique(stats.file.lg$lg),
																	 n.sites = n.sites,
																	 num.outliers = num.outliers,
																	 nnd.mean.null = null.mean, 
																	 nnd.sd.null = sd(null.mean.nnds, na.rm = TRUE),
																	 nnd.mean.emp = empirical.mean.nnd,
																	 nnd.emp.percentile = null.ecdf(empirical.mean.nnd),
																	 nnd.emp.zscore = (empirical.mean.nnd - null.mean)/sd(null.mean.nnds, na.rm = TRUE),
																	 nnd.emp.pvalue = two_side_p(null.mean.nnds, empirical.mean.nnd))
		}
	}
	return(do.call("rbind", nnd.stats))
}
