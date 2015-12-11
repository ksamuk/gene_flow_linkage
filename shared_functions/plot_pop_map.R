plot_pop_map <- function (region, pop.dat, manual = FALSE, adjust.lat = 0, adjust.lon = 0){
	
	region.dat.sub <- region.dat %>%
		filter(Region == region)
	
	if (manual == FALSE){
		adjust.lat <- region.dat.sub$adjust.lat
		adjust.lon <- region.dat.sub$adjust.lon
	}
	
	lat <- c(region.dat.sub$lat.min-adjust.lat, region.dat.sub$lat.max+adjust.lat) 
	lon <- c(region.dat.sub$lon.min-adjust.lon, region.dat.sub$lon.max+adjust.lon)
	# the map itself
	map("worldHires", ylim = lat, 
			xlim = lon,
			fill = TRUE,
			col = "grey",
			resolution = 0, mar = c(1,1,1,1))
	box()
	
	# population points
	
	region.dat <- pop.dat %>% 
		filter(Region == region) %>%
		filter(N > 1)
	
	region.coords <- list(x = region.dat$longitude, y = region.dat$latitude)
	points(region.coords, pch = 20, cex = 3, col = as.factor(region.dat$ecotype))
	map.centre <- c(mean(region.dat.sub$lat.min, region.dat.sub$lat.max), mean(region.dat.sub$lon.min, region.dat.sub$lon.max))
	map.scale(ratio = FALSE)
	#map.scale(map.centre[2], map.centre[1], units = "Km")
	
}