calculate.null.nnd <- function(stats.file.lg, num.outliers){
	
	site.sample <- stats.file.lg %>%
		filter(!is.na(gen.pos)) %>%
		select(gen.pos) %>%
		sample_n(num.outliers) %>%
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
	return(mean(nn.dist, na.rm = TRUE))
}