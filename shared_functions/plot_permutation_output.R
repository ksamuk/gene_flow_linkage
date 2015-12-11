plot_permutation_output <- function(permutation_output, stat, pal = pal, theme_all = theme_all){
	names(permutation_output$permuted.means.df)[1] <- "group"
	names(permutation_output$ecdf.df)[1] <- "group"
	
	plot.permute <- ggplot(data = permutation_output$permuted.means.df, aes_string(x = stat, fill = "group")) +
		geom_histogram() +
		geom_segment(data = permutation_output$ecdf.df ,aes(x = obs, xend = obs, y = 1,yend = 0, show_guide = F), 
								 size = 1, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed"))+
		scale_fill_manual(values = pal)+
		theme_all+
		facet_grid(~group, scales = "free")
	
}

