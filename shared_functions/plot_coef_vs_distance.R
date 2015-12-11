plot_coef_vs_distance <- function(all.df, coeff_type, coeff_name, size_all, alpha_all, theme_all){
	
	all.df <- all.df %>%
		mutate_(recomb_rate = interp(coeff_type))
	
	isolation_sum_plot <- all.df %>%
		filter(n_windows_fst > 100) %>%
		filter(pop1 != "lq") %>%
		mutate(comp_isolation = isolation.sum + least.cost.distance) %>%
		ggplot(aes(x =  sqrt(isolation.sum+1), y = log(recomb_rate_fst+3), color = group2, label = comparison)) +
		geom_point(size = 3, alpha = 1)+
		scale_color_manual(values = c(pal,1))+
		geom_smooth(method = "lm", size = 2)+
		#geom_smooth(data = all.df, aes(x = log(isolation.sum + least.cost.distance + 1), y = log(recomb_rate+3), color = "zblack"), size = 2)+
		facet_wrap(~ ecology)+
		theme_all+
		xlab("Non-oceanic distance (km)")+
		ylab(coeff_name)

all.df$euc.distance %>% sqrt %>% hist
	
	
	gc_distance_plot <- all.df %>%
		ggplot(aes(x = euc.distance, y = recomb_rate, color = ecology, label = comparison))+
		geom_point(size = size_all, alpha = alpha_all)+
		scale_color_manual(values = pal[3:4])+
		geom_smooth(method = "lm", size = size_all)+
		theme_all+
		facet_wrap(~ ecology)+
		xlab("Great circle distance (km)")+
		ylab("")
	
	lc_distance_plot <- all.df %>%
		ggplot(aes(x = least.cost.distance, y = recomb_rate, color = ecology, label = comparison)) +
		geom_point(size = size_all, alpha = alpha_all)+
		scale_color_manual(values = pal[3:4])+
		geom_smooth(method = "lm", size = size_all)+
		theme_all+
		facet_wrap(~ ecology)+
		xlab("Least cost distance (km)")+
		ylab(coeff_name)
	
	fst_plot <- all.df %>%
		ggplot(aes(x = fst, y = recomb_rate, color = ecology, label = comparison)) +
		geom_point(size = size_all, alpha = alpha_all)+
		#geom_text(size = 3)+
		scale_color_manual(values = pal[3:4])+
		geom_smooth(method = "lm", size = size_all)+
		theme_all +
		facet_wrap(~ ecology)+
		xlab("Average FST")+
		ylab("")
	
	plot_grid(isolation_sum_plot, gc_distance_plot, lc_distance_plot, fst_plot, align = "hv")
	
}