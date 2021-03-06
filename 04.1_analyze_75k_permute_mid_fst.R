
rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("ggplot2")
library("dplyr")
library("wesanderson")
library("grid")
library("gridExtra")
library("ggthemes")
library("lazyeval")
library("Hmisc")
library("cowplot")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select

################################################################################
# plotting parameters + global variables
################################################################################

theme_all <- theme_classic(base_size = 12)+
	theme(strip.text.x = element_blank(), 
				strip.background = element_blank(), 
				axis.title.x = element_blank(), 
				axis.title.y = element_blank(),
				#plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
				plot.margin = unit(c(0,0,0,0),"cm"),
				legend.position = "none",
				axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

#pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
pal <- c("#F21A00", "#3B9AB2" , "#E7C11A", "#9BBD95")

#number of permutations for tests
n_permutations <- 10000
options(scipen = -1)

################################################################################
# initialize coefficient data files
################################################################################

coeff.dat <- initialize_coeff_dat_files()

# get genome-wide fst values
fst.df <- read.table("meta_data/fst_df.txt", header = TRUE, stringsAsFactors = FALSE)

coeff.dat <- left_join(coeff.dat, fst.df)
	
################################################################################
# perform permutations + plot results
################################################################################

# remove comparisons that have fewer than 100 total shared windows (reduces noise)
# only retain comparisons with intermediate fsts (unique to this analysis)
coeff.dat <- coeff.dat %>%
	filter(n_windows_fst > 100) %>%
	filter(fst >= 0.185, fst <= 0.675)
	
# build a list of combinations of groups (strict, relaxed) 
# and fit types (fst, fst_dxy, dxy) to permute over (and plot / save pvalues)

stat <- c("recomb_rate_fst", "recomb_rate_dxy", "recomb_rate_hs")
stat <- rep(stat, each = 1)
group_type <- c("group2.new")
group_type <- rep(group_type, 3)
combo.df <- data.frame(stat, group_type, stringsAsFactors = FALSE)

# interate over combo.df (permute, extract p-values, plot)

# initialize results lists
pvals <- list()
plots <- list()
permutation_output <- list()

for (i in 1:nrow(combo.df)){
	
	permutation_output[[i]] <- run_cluster_permutations(coeff.dat, combo.df$stat[i], combo.df$group_type[i], n_permutations)
	pvals[[i]] <- save_pvals(permutation_output[[i]], combo.df$stat[i], combo.df$group_type[i])

}

for (i in 1:nrow(combo.df)){
	
	plots[[i]] <- plot_permutation_output(permutation_output[[i]], combo.df$stat[i], pal = pal, theme_all = theme_all)
	
}

################################################################################
# write pvalue table to file
################################################################################
pval.out <- pvals %>% bind_rows 
pval.out$group <- gsub("\\n"," ",pval.out$group)
pval.out$group_type <- gsub("group2.new", "relaxed", pval.out$group_type)

pval.out <- pval.out %>%
	mutate(effect.size = as.numeric((observed.means - permuted.means)))
	

write.table(pval.out, file = "meta_data/recomb_permutation_pvalues_mid_fst.txt", row.names = FALSE, quote = FALSE)

################################################################################
# plot figure S2
################################################################################

# note this plots the "raw" figure, annotations seen in paper were made manually

# similar to figure 1C 
fst_relaxed <- coeff.dat %>% 
	filter(!is.na(recomb_rate_fst)) %>%
	filter(n_windows_fst > 100) %>%
	plot_dot_line_plot(., group = "group2.new", stat = "recomb_rate_fst", label = "", 
										 pal = pal, y_lab = expression('F'["ST"]*" Outlier vs. Recombination Coefficient"), theme_all = NULL, 
										 point_size = 1, line_size = 2)

dxy_relaxed <- coeff.dat %>% 
	filter(n_windows_dxy > 100) %>% 
	filter(!is.na(recomb_rate_dxy)) %>%
	plot_dot_line_plot(., group = "group2.new", stat = "recomb_rate_dxy", label = "", 
										 pal = pal, y_lab = "", theme_all = NULL, 
										 point_size = 1, line_size = 2)

hs_relaxed <- coeff.dat %>% 
	filter(!is.na(recomb_rate_hs)) %>%
	filter(n_windows_fst > 100) %>%
	plot_dot_line_plot(., group = "group2.new", stat = "recomb_rate_hs", label = "", 
										 pal = pal, y_lab = "", theme_all = NULL, 
										 point_size = 1, line_size = 2)

labels <- c("fst", "dxy", "hs")
figS2_dot <- plot_grid(fst_relaxed, dxy_relaxed, hs_relaxed, 
									 ncol = 1, labels = labels, align = "hv")

save_plot(figS2_dot, filename = "figures/raw/FigureS2_mid_fst_dotplots_raw.pdf", base_height = 8.5, base_width = 4.25)

figS2_hist <- plot_grid(plots[[1]], plots[[2]], plots[[3]], 
									 ncol = 1, labels = labels, align = "hv")

save_plot(figS2_hist, filename = "figures/raw/FigureS2_mid_fst_hist_raw.pdf", base_height = 8.5, base_width = 4.25)

plot_grid(fst_relaxed, plots[[1]], dxy_relaxed, plots[[2]], hs_relaxed, plots[[3]], 
					ncol = 2, labels = labels, align = "hv")

