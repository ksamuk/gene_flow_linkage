################################################################################
# This script performs the following functions:
# 1. Plots all the raw figures for the manuscript
# 3. Performs permutations on fitted model ouput
# 2. Creates all the supplementary data files
################################################################################

library("dplyr")

# this generates warnings; they are safe to ignore (mostly ggplot reporting NAs being removed)
list.files(pattern = "plot", full.names = TRUE) %>% 
	grep("RUN", ., invert = TRUE, value = TRUE) %>% 
	sapply(source)

list.files(pattern = "plot", full.names = TRUE) %>% 
	grep("RUN", ., invert = TRUE, value = TRUE) %>% 
	sapply(source)
