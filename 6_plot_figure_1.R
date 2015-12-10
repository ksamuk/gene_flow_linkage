# figure 1 
rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("ggplot2")
library("readr")
library("dplyr")
library("wesanderson")
library("Hmisc")
library("ggthemes")
library("gridExtra")
library("cowplot")
library("scales")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

################################################################################
# plot color scheme/theme
################################################################################

#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

theme_all <- theme_hc(base_size = 16) +
	theme(axis.title.y = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				plot.margin = unit(c(0,0,0,0),"cm"))

################################################################################
# load raw data
################################################################################

coeff.dat <- initialize_coeff_dat_files()

stats_df <- read_delim(file = "analysis_ready/75k_stats_combined.txt", delim = " ")
stats_df <- add_region_data(stats_df)

stats_df <- stats_df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = "."))

rep.comparisons <- c("cp.marine.sk.marine", 
										 "constance.stream.joes.stream", 
										 "mariager.marine.joes.lake", 
										 "pri.limnetic.pri.benthic")

stats_df <- stats_df %>% 
	filter(comparison %in% rep.comparisons)

################################################################################
# Figure 1 A: Representitive FST distributions
################################################################################

stats_df_filt <- stats_df %>%
	filter(!is.na(fst)) %>%
	filter(!is.infinite(fst)) %>%
	filter(recomb_rate <= 25) %>% 
	filter(var.sites >= 2) %>%
	filter(ds <= 3) %>%
	filter(gene_count <= 15) %>%
	filter(lg!=19) %>%
	mutate(dxy_adj = ifelse(sites >= 500, dxy, NA)) %>% 
	group_by(comparison) %>%
	mutate(fst.outlier = is.outlier(fst))%>%
	mutate(dxy.outlier = is.outlier(dxy)) %>%
	mutate(both.outlier = fst.outlier & dxy.outlier)%>%
	mutate(hs = (hexp1+hexp2)/2) %>%
	ungroup

rep_plot <- stats_df %>%
	group_by(comparison) %>%
	filter(lg  == 4) %>%
	mutate(fst_outlier_value = ifelse(fst.outlier == TRUE, fst, NA)) %>%
	mutate(fst_non_outlier_value = ifelse(fst.outlier == FALSE, fst, NA)) %>%
	ungroup %>%
	ggplot(aes(x = midpos/1000000, y = fst))+
	geom_point(aes(x = midpos/1000000, y = fst_outlier_value, color = "zzzzred"))+
	geom_point(aes(x = midpos/1000000, y = fst_non_outlier_value, color = "zzzzgrey"))+
	stat_smooth(method = "loess", aes(color = group2), size = 1.5, se = FALSE)+
	stat_smooth(method = "loess", aes(color = "zzz", x = midpos/1000000, y = recomb_rate/10), size = 1.5, se = FALSE, linetype ="11111111")+
	theme_hc_border()+
	ylab(expression('F'["ST"]))+
	xlab("Chromosomal position (MB)")+
	ylim(0, 1)+
	theme(panel.grid = element_blank(),
				strip.background = element_blank(),
				legend.position = "none",
				panel.border = element_rect(color="grey", fill=NA),
				strip.text = element_blank())+
	facet_grid(group2~lg)+
	scale_color_manual(values = c(pal, 1, "grey" ,1))+
	scale_x_continuous(labels = comma)

ggsave(rep_plot, filename = "figures/raw/Figure1_A_raw.pdf", height = 8.5, width = 4.25)

################################################################################
# Figure 1 B: Averaged Fitted model coefficients 
################################################################################

# create functions based on the averaged the regression coefficients for each selection/gene flow group
avg_regression_functions <- coeff.dat %>% filter(n_windows_fst > 100) %>% create_average_regression_functions_fst
#avg_regression_functions <- coeff.dat %>% filter(n_windows_dxy > 200) %>% create_average_regression_functions_dxy

# plot the functions 

# an empty plot
avg_reg_plot <- ggplot()+
	scale_x_continuous(limits=c(0, 25), expand=c(0,0))+
	scale_y_continuous(limits=c(0, 0.175), expand = c(0,0)) +
	xlab("Recombination rate (cM/MB)") +
	ylab(expression('F'["ST"]*" Outlier Probability"))

# add in the avg reg function plots 
for(i in 1:length(avg_regression_functions)){
	avg_reg_plot <- avg_reg_plot + 
		stat_function(aes(y = 0), fun = avg_regression_functions[[i]], colour = pal[i], size = 3)
}
	

################################################################################
# Figure 1 C: Jitter plot of model coefficients 
################################################################################

fst_relaxed <- coeff.dat %>% 
	filter(!is.na(recomb_rate_fst)) %>%
	filter(n_windows_fst > 100) %>%
	plot_dot_line_plot(., group = "group2.new", stat = "recomb_rate_fst", label = "", 
										 pal = pal, y_lab = expression('F'["ST"]*" Outlier vs. Recombination Coefficient"), theme_all = NULL, 
										 point_size = 1, line_size = 2)

################################################################################
# Combine figures
################################################################################

figure_1_plot <- plot_grid(rep_plot, avg_reg_plot, NULL,fst_relaxed, align = "hv", ncol = 2, rel_heights = c(1,1))

save_plot("figures/raw/Figure1_B_C_raw.pdf", figure_1_plot, base_height = 8.5, base_width = 8.5)

# Figure1_A_raw.pdf & Figure1_B_C_raw.pdf are then combined/cleaned up manually in illustrator 

################################################################################
# Supplemental figures
################################################################################

# # Supp matt
# dxy_relaxed <- coeff.dat %>% 
# 	filter(n_windows_dxy > 100) %>% 
# 	filter(!is.na(recomb_rate_dxy)) %>%
# 	plot_dot_line_plot(., group = "group2.new", stat = "recomb_rate_dxy", label = "", 
# 										 pal = pal, y_lab = "", theme_all = NULL, 
# 										 point_size = 1, line_size = 2)
# 
# hs_relaxed <- coeff.dat %>% 
# 	filter(!is.na(recomb_rate_hs)) %>%
# 	filter(n_windows_fst > 100) %>%
# 	plot_dot_line_plot(., group = "group2", stat = "recomb_rate_hs", label = "", 
# 										 pal = pal, y_lab = "", theme_all = NULL, 
# 										 point_size = 1, line_size = 2)


