## permutation test to assess significance of clustering data

rm(list = ls())

################################################################################
# Libraries and shared functions
################################################################################

library("ggplot2")
library("dplyr")
library("wesanderson")
library("grid")
library("gridExtra")
library("ggthemes")
library("lazyeval")
library("Hmisc")

list.files("shared_functions", full.names = TRUE) %>% sapply(source)

select <- dplyr::select

################################################################################
# Read in and format input files
################################################################################

cluster.df <- initialize_clustering_output()

# filter out linakge groups with super low data
# clustering metrics behave strangely when number of sites / number of outliers is very low
cluster.df <- cluster.df %>%
	filter(n.sites > 30) %>%
	filter(num.outliers > 5)

pvals <- list()
plots <- list()

theme_all <- theme_classic(base_size = 12)+
		theme(strip.text.x = element_blank(), 
				strip.background = element_blank(), 
				axis.title.x = element_blank(), 
				axis.title.y = element_blank(),
				#plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
				plot.margin = unit(c(0,0,0,0),"cm"),
				legend.position = "none"
				#axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
				)

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

n_permutations <- 10000

################################################################################
# Permutation test : nnd.diff, relaxed
################################################################################

stat <- "nnd.diff.sd"
group_type <- "group2.new" # relaxed

permutation_output <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group2.new, comparison) %>%
	summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	filter(!is.na(nnd.diff.sd)) %>%
	filter(nnd.diff.sd > -4)%>% # remove single crazy low value (some kind of error)
	run_cluster_permutations(., stat, group_type, n_permutations, collapsed = TRUE)

pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# Permutation test : dispersion, relaxed
################################################################################

stat <- "disp.out"
group_type <- "group2.new" # relaxed

permutation_output <- run_cluster_permutations(cluster.df, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# make unified plot
################################################################################

# theme settings
point_size <- 1
line_size <- 3
size <- 16
theme.all <- theme_hc + theme(legend.position="none", 
															axis.title.x = element_blank(), 
															axis.title.y = element_text(vjust=1.5),
															axis.text.x = element_blank(),
															axis.line.x = element_blank(),
															axis.ticks.x = element_blank(),
															strip.text = element_blank(), 
															strip.background = element_blank(),
															legend.title = element_blank())

# coeff of disp
coeff.dispersion <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group2.new, comparison) %>%
	summarise(disp.out = mean(disp.out)) %>%
	plot_dot_line_plot(data = ., group = "group2.new", 
										 stat = "disp.out", label = "", 
										 pal = pal, y_lab = "Outlier dispersion coefficient",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)

# nnd.diff
nnd.diff <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group2.new, comparison) %>%
	summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	#summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	filter(!is.na(nnd.diff.sd)) %>%
	filter(nnd.diff.sd > -4)%>% # filters out one particular impossible value
	plot_dot_line_plot(data = ., group = "group2.new", 
										 stat = "nnd.diff.sd", label = "", 
										 pal = pal, y_lab = "Expected NND - Outlier NND (cM)",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)


labels <- c("nnd.diff_relaxed", "disp_relaxed")
figS3<- plot_grid(nnd.diff, plots[[1]], coeff.dispersion, plots[[2]], ncol = 2, labels = labels, align = "hv")
save_plot(figS3, filename = "figures/raw/FigureS3_raw.pdf", base_height = 8.5, base_width = 8.5)

################################################################################
# format and write pvalue table to file
################################################################################

pval.out <- pvals %>% bind_rows
pval.out$group <- gsub("\\n"," ",pval.out$group)
pval.out[pval.out$group_type=="group2.new",]$group_type <- "relaxed"
pval.out[pval.out$group_type=="group.new",]$group_type <- "strict"
pval.out <- pval.out %>%
	mutate(effect.size = observed.means - permuted.means)

write.table(pval.out, file = "meta_data/clustering_permutation_pvalues_fst.txt", row.names = FALSE, quote = FALSE)


