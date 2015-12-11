plot_dot_line_plot <- function(data, group, stat, label = "", pal = pal, y_lab = "", theme_all = NULL, point_size = 1, line_size = 2){
	
	if (length(theme_all) <= 0){
		theme_all <- theme_hc(base_size = 16)
			#+theme(plot.margin = unit(c(1,1,1,1),"cm"))
	}
	
	plot <- data %>%
		ggplot(aes_string(x = group, y = stat))+
		geom_point(position = position_jitter(width = 0.3), size = point_size)+
		geom_errorbar(stat = "hline", 
									yintercept = "mean",
									width = 0.8, size = line_size, 
									aes(ymax = ..y.. , ymin = ..y..),
									color = c(pal[1], pal[2], pal[3], pal[4]))+
		#geom_text(mapping = NULL, aes_string(label = label, hjust = 0, vjust = 0), color = "black" ,alpha = 0.25, size = 6)+
		theme_all+
		ylab(y_lab)+
		xlab(NULL)
	
	return(plot)
	
}