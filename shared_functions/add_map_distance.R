add_map_distance <- function(snp.file){
	
	# read in snp file
	
	snp.file$gen.pos <- NA
	snp.file$pos <- as.numeric(snp.file$pos)
	
	snp.file <- snp.file %>%
		arrange(lg, pos)
	
	# read in genetic map data
	map.file <- list.files(file.path("ev_prep_scripts"), pattern="roesti_recomb_estimates.txt", full.names = TRUE)
	map.file <- read.table(map.file, header = TRUE)
	names(map.file)[1] <- "lg"
	
	
	# calculate genetic distance
	
	#place holder vector
	
	snp.map.position <- c()
	
	for (i in unique(snp.file$lg)){
		snp.file.lg <- snp.file %>%
			filter(lg == i)
		
		map.file.lg <- map.file %>%
			filter(lg == i)
		
		snp.map.position.lg <- snp.file.lg$gen.pos
		
		for (k in 1:length(snp.file.lg$pos)){
			
			snp.pos <- snp.file.lg$pos[k]
			
			#snp.pos <- 1190613 - 12500
			
			# find matching window(s) for snp
			# because of the reversed sections in the map file (Roesti et al. 2013)
			# this can be 1 (normal condition),2 (position is only found in breakpoints of reverse section), 
			# or 3 (found in breakpoints and middle of reversed section)
			windows <- bind_rows(map.file.lg[map.file.lg$pos1 < snp.pos & map.file.lg$pos2 >= snp.pos,],
									 map.file.lg[map.file.lg$pos1 > snp.pos & map.file.lg$pos2 <= snp.pos,])
			
			snp.pos.window <- windows
			
			if (nrow(windows) == 3){
				#print(paste(i, snp.pos, nrow(windows)))
				# marker is in a reversed section
				window.diff <- (windows$pos1 - windows$pos2) %>% abs
				# choose window with smallest physical distance (i.e. not the breakpoints)
				snp.pos.window <- windows[which(window.diff == min(window.diff)),]
			} else if (nrow(windows) == 2){
				#print(paste(i, snp.pos, nrow(windows)))
				# marker only found in reversed section breakpoints, set as blank
				snp.pos.window <- c()
			}
		
			if (length(snp.pos.window$lg) == 1){
				
				#ensure markers are in order
				if (snp.pos.window$pos1 > snp.pos.window$pos2){
					snp.pos.window <- with(snp.pos.window, data.frame(lg = lg, pos1 = pos2, pos2 = pos1, gen1 = gen2, gen2 = gen1))
				}
				
				# how far (proportionally) snp is along window in bp
				snp.window.position <- 1 - (snp.pos.window$pos2 - snp.pos) / (snp.pos.window$pos2 - snp.pos.window$pos1)
				
				# convert to genetic distance
				snp.map.position.lg[k] <- (snp.window.position * (snp.pos.window$gen2 - snp.pos.window$gen1)) + snp.pos.window$gen1
			} else{
				snp.map.position.lg[k] <- NA
			}
			
			
		}
		
		snp.map.position<-c(snp.map.position, snp.map.position.lg)
		
	}
	
	snp.file$gen.pos <- snp.map.position
	return(snp.file)
}