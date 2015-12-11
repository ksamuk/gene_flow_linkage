library(ggplot2)

visualize_75k_stats_file <- function (stats_file, stat, plot_type){
	
	stats_file <- read.table(stats_file, header = TRUE, 
													 stringsAsFactors = FALSE)
	
	stats_file$Chr <- chrom.to.num(stats_file$Chr)
	stats_file <- stats_file %>%
		mutate(fst.outlier = is.outlier(Fst)) %>%
		mutate(dxy.outlier = is.outlier(Dxy))
	
	stats_file %>%
		ggplot(aes(x = EndPos, y = Fst, color = fst.outlier))+
		geom_point()+
		facet_wrap(~Chr)
	
	stats_file %>%
		ggplot(aes(x = EndPos, y = Dxy, color = dxy.outlier))+
		geom_point()+
		facet_wrap(~Chr)
		
}	