# plot figure 2

rm(list=ls())

################################################################################
# Libraries and shared functions
################################################################################

library("ggplot2")
library("dplyr")
library("wesanderson")
library("grid")
library("gridExtra")
library("ggthemes")
library("cowplot")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible()

select <- dplyr::select
#pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
pal <- c("#F21A00", "#3B9AB2" , "#E7C11A", "#9BBD95")

################################################################################
# Read in and format input files
################################################################################

# the clustering file
cluster.df <- initialize_clustering_output()

# filter out linakge groups with super low data
cluster.df <- cluster.df %>%
	filter(n.sites > 30) %>%
	filter(num.outliers > 5)

################################################################################
# Figure 2 (clustering data, relaxed)
################################################################################

# theme settings
point_size <- 1
line_size <- 1
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

cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group2.new, comparison) %>%
	summarise(disp.out = mean(disp.out)) %>%
	ungroup %>%
	ggplot(aes_string(x = "group2.new", y = "disp.out"))+
	geom_point(position = position_jitter(width = 0.3), size = point_size, color = "grey50", alpha = 0.5)+
	stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
							 geom = "crossbar", width = 0.5, size =1, color = c(pal[1], pal[2], pal[3], pal[4]))


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

fig2 <- plot_grid(coeff.dispersion, nnd.diff)
save_plot(fig2, filename = "figures/raw/Figure2_raw.pdf", base_height = 4.5, base_width = 8.5)
