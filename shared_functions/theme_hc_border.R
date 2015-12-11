theme_hc_border <- function (base_size = 12, base_family = "sans", bgcolor = "default") 
{
	bgcol <- ggthemes_data$hc$bg[bgcolor]
	ret <- theme(text = element_text(size = base_size, family = base_family), 
							 title = element_text(hjust = 0.5), axis.title.x = element_text(hjust = 0.5), 
							 axis.title.y = element_text(hjust = 0.5), panel.grid.major.y = element_line(color = "gray"), 
							 panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), 
							 panel.grid.minor.x = element_blank(), 
							 panel.background = element_blank(), legend.position = "bottom", 
							 legend.key = element_rect(fill = "#FFFFFF00"))
	ret
}