rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("mapdata")
library("maps")
library("maptools")
library("rworldmap")
library("dplyr")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

################################################################################
# Plot figure S1
################################################################################

# initialize the world map data

map("worldHires",  xlim = c(-130, -100), ylim = c(50,60), resolution = 0)

# input files
pop.dat <- read.csv(file = "meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
pop.dat <- pop.dat %>%
	select(pop, ecotype, latitude, longitude, Region, N) %>%
	filter(N > 1) %>%
	filter(Region != "Russia")

# some tweaks for 'zoomed in' maps

region.dat <-  pop.dat %>%
	group_by(Region) %>%
	summarise(lat.min = min(latitude), lat.max = max(latitude), 
						lon.min = min(longitude), lon.max = max(longitude))

region.dat$adjust.lat <- c(2.05, 1.25, 2.55, 17, 1.9, 2, 5)
region.dat$adjust.lon <- c(6, 1, 4.5, 25, 2.5, 2, 8)

# plot individual maps to file
pdf(file = "figures/raw/figureS1_raw.pdf", width = 8.5, height = 8.5)
par(mfrow = c(3,3), mar = c(0,0,0,0))
lapply(region.dat$Region, plot_pop_map, pop.dat = pop.dat)
dev.off()
