############################################################
# Plot Figure 3
# KS Nov 2015
############################################################

rm(list =ls())

#########################
# LIBRARIES & FUNCTIONS
#########################

library("dplyr")
library("tidyr")
library("broom")
library("ggplot2")
library("ggthemes")
library("Hmisc")
library("cowplot")
library("lazyeval")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

#########################
# INPUT FILES
#########################

# df containing mean fst values for all populations
fst.df <- read.table("meta_data/fst_df.txt", header = TRUE, stringsAsFactors = FALSE)

# dfs containing linear model coefficients for outliers vs. recomb rate
fst.model.fits <- initialize_coeff_dat_files()

# df containing distances (least cost, great circle) between populations
# also measures of isolation from ocean (dist to coast, i.e. how inland, for each population)
dist.df <- read.table("meta_data/pop_geo_distances.txt", header = TRUE, stringsAsFactors = FALSE)

#########################
# PRE PROCESSING
#########################

# bind the three data files together
all.df <- cbind(fst.model.fits, dist.df[,7:10], fst.df$fst)
names(all.df)[length(all.df)] <- "fst"

# add in combined measures of allopatry (simple transformations of distance data)
all.df <- all.df %>% 
	mutate(isolation.sum = dist.to.coast1 + dist.to.coast2)%>%
	mutate(comparison = paste0(pop1, pop2))%>%
	rowwise() %>%
	mutate(isolation.max = max(dist.to.coast1, dist.to.coast2))

all.df$ecology <- gsub("D", "Divergent", all.df$ecology)
all.df$ecology <- gsub("S", "Parallel", all.df$ecology)

##################################
# Size / theme for figs s7 - s9
##################################

size_all <- 2
alpha_all <- 0.2

theme_all <- theme_hc(base_size = 16) + 
	theme(legend.position = "none",
				strip.background = element_blank()) 

##################################
# Figure 3
##################################

# recomb rate vs. least cost distance

dist_fst <- all.df %>%
	filter(n_windows_fst > 200) %>%
	ggplot(aes(x =  jitter(sqrt(euc.distance +1), 2000), y = recomb_rate_fst, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 1)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, color = "black", se = TRUE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab(expression(sqrt("Great circle distance (km)")))+
	ylab("Recombination rate vs. \nFST outlier coeffficient")

fst_fst <- all.df %>%
	filter(n_windows_fst > 200) %>%
	ggplot(aes(x =  fst, y = recomb_rate_fst, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 1)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, se = TRUE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab("Genome-wide average FST")+
	ylab("Recombination rate vs. \nFST outlier coeffficient")

fig3 <- plot_grid(dist_fst, fst_fst, ncol = 1)

save_plot(fig3, file = "figures/raw/Figure3_raw.pdf", base_height = 8.5, base_width = 8.5)

##################################
# Supporting linear models
##################################

# fst: distance 

fit_dist_bias <- all.df %>%
	filter(n_windows_fst > 200) %>%
	mutate(sqrt_euc = sqrt(euc.distance+1)) %>%
	lm(data = ., recomb_rate_fst~sqrt_euc * ecology)

anova(fit_dist_bias)
summary(fit_dist_bias)

# fst: fst*ecology

fit_ecol_bias <- all.df %>%
	filter(n_windows_fst > 200) %>%
	mutate(sqrt_euc = sqrt(euc.distance+1)) %>%
	lm(data = ., recomb_rate_fst~fst*group2)

anova(fit_ecol_bias)
summary(fit_ecol_bias)

# dxy: distance 

fit_dist_bias_dxy <- all.df %>%
	filter(n_windows_dxy > 200) %>%
	mutate(sqrt_euc = sqrt(euc.distance+1)) %>%
	lm(data = ., recomb_rate_dxy ~ sqrt_euc)

anova(fit_dist_bias_dxy)
summary(fit_dist_bias_dxy)

# dxy: dxy*ecology

fit_ecol_bias_dxy <- all.df %>%
	filter(n_windows_dxy > 200) %>%
	mutate(sqrt_euc = sqrt(euc.distance+1)) %>%
	lm(data = ., recomb_rate_dxy~fst*group2)

anova(fit_ecol_bias_dxy)
summary(fit_ecol_bias_dxy)

##################################
# Figure S4
##################################

dist_fst <- all.df %>%
	filter(n_windows_dxy > 200) %>%
	ggplot(aes(x =  jitter(sqrt(euc.distance +1), 2000), y = recomb_rate_dxy, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 1)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, color = "black", se = TRUE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab(expression(sqrt("Great circle distance (km)")))+
	ylab("Recombination rate vs. \nDxy outlier coeffficient")

fst_fst <- all.df %>%
	filter(n_windows_dxy > 200) %>%
	ggplot(aes(x =  fst, y = recomb_rate_dxy, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 1)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, se = TRUE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab("Genome-wide average FST")+
	ylab("Recombination rate vs. \nDxy outlier coeffficient")


figs4 <- plot_grid(dist_fst, fst_fst, ncol = 1)

save_plot(fig3, file = "figures/raw/FigureS4_raw.pdf", base_height = 8.5, base_width = 8.5)
