############################################################
# Plot Figure S4 & S5
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
# Size / theme for figs s4 & s5
##################################

size_all <- 2
alpha_all <- 0.2

theme_all <- theme_hc(base_size = 16) + 
	theme(legend.position = "none",
				strip.background = element_blank()) 

##################################
# Figure S4
##################################

# recomb rate vs. least cost distance

dist_fst <- all.df %>%
	filter(n_windows_fst > 200) %>%
	ggplot(aes(x =  jitter(sqrt(euc.distance +1), 2000), y = recomb_rate_fst, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 0.5)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, color = "black", se = FALSE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab(expression(sqrt("Great circle distance (km)")))+
	ylab("Recombination rate vs. \nFST outlier coeffficient")

fst_fst <- all.df %>%
	filter(n_windows_fst > 200) %>%
	ggplot(aes(x =  fst, y = recomb_rate_fst, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 0.5)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, se = FALSE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab("Genome-wide average FST")+
	ylab("Recombination rate vs. \nFST outlier coeffficient")

figs4 <- plot_grid(dist_fst, fst_fst, ncol = 1)

save_plot(figs4, file = "figures/raw/FigureS4_raw.pdf", base_height = 8.5, base_width = 8.5)

##################################
# Figure S5
##################################

dist_dxy <- all.df %>%
	filter(n_windows_dxy > 200) %>%
	ggplot(aes(x =  jitter(sqrt(euc.distance +1), 2000), y = recomb_rate_dxy, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 0.5)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, color = "black", se = FALSE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab(expression(sqrt("Great circle distance (km)")))+
	ylab("Recombination rate vs. \nDxy outlier coeffficient")

fst_dxy <- all.df %>%
	filter(n_windows_dxy > 200) %>%
	ggplot(aes(x =  fst, y = recomb_rate_dxy, color = group2, label = comparison)) +
	geom_point(size = 2, alpha = 0.5)+
	scale_color_manual(values = c(pal,1))+
	geom_smooth(method = "lm", size = 2, se = FALSE)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab("Genome-wide average FST")+
	ylab("Recombination rate vs. \nDxy outlier coeffficient")


figs5 <- plot_grid(dist_dxy, fst_dxy, ncol = 1)

save_plot(figs5, file = "figures/raw/FigureS5_raw.pdf", base_height = 8.5, base_width = 8.5)

##################################
# Supporting linear models (NB: not used, replaced by permutation tests)
##################################

# DISTANCE VS. RECOMBINATION BIAS

dist_df <- all.df %>%
	select(comparison, geography2, ecology, group2, recomb_rate_fst, euc.distance)

# confirming groups are correct
dist_df %>%
	ggplot(aes(x = euc.distance, y = recomb_rate_fst, color = ecology))+
	geom_point() +
	geom_smooth(method = "lm")

# permutation to test for difference in slopes
# computes the expected difference between selective regimes in the slope of the relationship between 
# recombination bias and geographic distance

permute_slope_distance <- function(data_df){
	
	perm_df <- dist_df
	perm_df$ecology <- sample(perm_df$ecology)
	
	perm_df <- perm_df %>%
		group_by(ecology) %>%
		do(slope = lm(.$recomb_rate_fst ~ .$euc.distance) %>% coefficients %>% .[2] %>% unlist) %>%
		data.frame
	
	perm_df$slope <- unlist(perm_df$slope)
	
	
	perm_df$slope[1] - perm_df$slope[2]
	
}

# perform the permutation
perm_reps <- replicate(10000, permute_slope_distance(dist_df), simplify = FALSE)
perm_reps <- unlist(perm_reps)

# compute the empirical difference in means
emp_val <- dist_df %>%
	group_by(ecology) %>%
	do(slope = lm(.$recomb_rate_fst ~ .$euc.distance) %>% coefficients %>% .[2] %>% unlist) %>%
	data.frame

emp_val$slope <- unlist(emp_val$slope)
emp_val <- emp_val$slope[1] - emp_val$slope[2]

# compute a two-sided pvalue by determining the location of the empirical value in the distribution of 
# permuted values
two_side_p(perm_reps, emp_val)

# 10000 permutations:
#[1] 0.00079992

# GENOME-WIDE FST VS. RECOMBINATION BIAS

fst_df <- all.df %>%
	select(comparison, geography2, ecology, group2, recomb_rate_fst, fst)

# confirming groups are correct
fst_df %>%
	ggplot(aes(x = fst, y = recomb_rate_fst, color = ecology))+
	geom_point() +
	geom_smooth()

# permutation to test for difference in slopes
# computes the expected difference between selective regimes in the slope of the relationship between 
# recombination bias and geographic distance

permute_slope_distance <- function(data_df){
	
	perm_df <- dist_df
	perm_df$ecology <- sample(perm_df$ecology)
	
	perm_df <- perm_df %>%
		group_by(ecology) %>%
		do(slope = lm(.$recomb_rate_fst ~ .$euc.distance) %>% coefficients %>% .[2] %>% unlist) %>%
		data.frame
	
	perm_df$slope <- unlist(perm_df$slope)
	
	
	perm_df$slope[1] - perm_df$slope[2]
	
}

# perform the permutation
perm_reps <- replicate(10000, permute_slope_distance(dist_df), simplify = FALSE)
perm_reps <- unlist(perm_reps)

# compute the empirical difference in means
emp_val <- dist_df %>%
	group_by(ecology) %>%
	do(slope = lm(.$recomb_rate_fst ~ .$euc.distance) %>% coefficients %>% .[2] %>% unlist) %>%
	data.frame

emp_val$slope <- unlist(emp_val$slope)
emp_val <- emp_val$slope[1] - emp_val$slope[2]

# compute a two-sided pvalue by determining the location of the empirical value in the distribution of 
# permuted values
two_side_p(perm_reps, emp_val)

# 10000 permutations:
#[1] 0.00079992
