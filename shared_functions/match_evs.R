match_evs <- function(stats.file.name, linear_model_function = linear_model_function){
	
	print(paste0("processing ", stats.file.name,"..."))
	#initialize data to be matched 
	stats.file <- read.table(stats.file.name, header = TRUE)
	if (length(stats.file[,1]) < 1){
		return(NA)
	}
	names(stats.file) <- c("lg","pos1","pos2","midpos","sites","var.sites","dxy","fst","hexp1","hexp2")
	stats.file$lg <- chrom.to.num(stats.file$lg)
	stats.file$lg <- as.numeric(stats.file$lg)
	matched.all <- data.frame(stats.file)
	###looping through all the ev files
	for (i in 1:length(ev.files)){
		
		#read ev file
		ev <- read.table(file = ev.files[i], header=TRUE)
		
		#prep empty ev dataframe
		matched.evs <- data.frame()
		length(matched.evs) <- 4
		names(matched.evs) <- c("lg", "pos1", "pos2", "ev")
		
		#sanitize file of NAs and duped rows
		ev <- ev[complete.cases(ev), ]
		ev <- arrange(ev, lg, pos1)
		ev <- unique(ev)
		
		###loop through lgs, matching evs as we god
		for (j in 1:max(stats.file$lg)){
			#for (j in 1:2){
			
			#subset ev and stat by lg
			ev.chr <- subset(ev,ev$lg == j)
			stat.chr <- subset(stats.file, stats.file$lg == j)
			stat.chr <- stat.chr[,1:3]
			
			#build IRanges objects for overlap finding
			ev.range <- IRanges(start = ev.chr$pos1, end = ev.chr$pos2)
			#minimal set. note that this means overlapping evs end up returning the same value, rather than say, their average
			#(this shouldn't matter for most types of evs)
			ev.range <- reduce(ev.range, min.gapwidth = 0L)
			stat.range <- IRanges(start = stat.chr$pos1, end = stat.chr$pos2)
			
			#find ovelaps amd build an "overlap df"
			overlap <- findOverlaps(stat.range, ev.range, select = "all")
			overlap.df <- data.frame(lg = stat.chr[queryHits(overlap),]$lg,
															 pos1 = stat.chr[queryHits(overlap),]$pos1,
															 pos2 = stat.chr[queryHits(overlap),]$pos2,
															 ev.start = start(ev.range[subjectHits(overlap)]),
															 ev.end = end(ev.range[subjectHits(overlap)]))
			
			#truncate overlaps that start before window
			overlap.df$ev.start[overlap.df$ev.start < overlap.df$pos1] <- overlap.df$pos1[overlap.df$ev.start < overlap.df$pos1]
			
			#truncate overlaps that extend past end of window
			overlap.df$ev.end[overlap.df$ev.end > overlap.df$pos2] <- overlap.df$pos2[overlap.df$ev.end > overlap.df$pos2]
			
			#calc widths
			overlap.df$width <- (overlap.df$ev.end - overlap.df$ev.start) + 1
			
			overlap.df$ev <- ev.chr[subjectHits(overlap), 4]
			overlap.df <- unique(overlap.df)
			
			#calculate the proportion of total overlap
			prop.overlap <- overlap.df %>%
				group_by(pos1) %>%
				summarise(total.overlap = sum(width))
			
			#add it back into the overlap dataframe
			prop.overlap <- data.frame(prop.overlap)
			overlap.df$total.overlap <- prop.overlap$total.overlap[match(overlap.df$pos1, prop.overlap$pos1)]
			
			#the weighted contribution this to the total (deals with multiple ev hits to a single window)
			overlap.df$ev.prop <- overlap.df$ev * (overlap.df$width / overlap.df$total.overlap)
			
			#calculate the weighted evs for each window
			weighted.evs <- overlap.df %>%
				group_by(pos1) %>%
				summarise(ev = sum(ev.prop))
			weighted.evs <- data.frame(weighted.evs)
			#add back into df
			overlap.df$ev.avg <- weighted.evs$ev[match(overlap.df$pos1, weighted.evs$pos1)]
			
			####DIVERGES FROM SNP BASED ANALYSIS
			#HACK TO ACCOMODATE WINDOWS re-read stats file
			stat.chr <- subset(stats.file, stats.file$lg==j)
			stat.chr$ev <- overlap.df$ev.avg[match(stat.chr$pos1, overlap.df$pos1)]
			matched.evs.chr <- stat.chr
			
			#rbind chromos
			matched.evs <- rbind(matched.evs, matched.evs.chr)
			
		}
		###end lg loop
		
		#attach real name of ev and cbind to stats file
		names(matched.evs)[11] <- sapply(strsplit(ev.files[i], split="/"), function(x)gsub(".txt","",x[length(x)]))
		matched.all <- suppressMessages(left_join(matched.all, matched.evs))
		
	}
	
	# filter and call outliers
	matched.all <- filter_data_call_outliers_function(matched.all)
	
	# if there is extremely low data, skip model fitting

		# fit linear model and catch failures to converge (e.g. due to insufficient data)
		
		#blank the warning/error messages
		linear_model_warning <- NULL
		
		# fit model and catch warning (if any)
		withCallingHandlers(model.out <- linear_model_function(matched.all), warning = function(w){
			linear_model_warning <<- w$message
			invokeRestart("muffleWarning")
		})
	
		# if model threw warning, print it and set coeffs as NA
		if(length(linear_model_warning) >= 1){
			print(linear_model_warning)
			#coeffs <- rep(list(NA), 4)
			#coeffs <- setNames(coeffs, c("intercept", "recomb_rate", "ds", "gene_count"))
		} else{
			linear_model_warning <- NA
		}
		coeffs <- as.list(model.out$coefficients)
		names(coeffs)[1] <- "intercept"
	
	# parse file name
	file.name.stripped <- sapply(strsplit(stats.file.name, split = "/"), function(x)gsub(".txt","",x[length(x)]))
	file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
	pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
	ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
	pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
	ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
	geography <- file.name.split[3]
	ecology <- file.name.split[4]
	n_windows_dxy <- sum(!is.na(matched.all$dxy.outlier))
	n_windows_fst <- sum(!is.na(matched.all$fst.outlier))
	
	row.out <- data.frame(pop1, ecotype1, pop2, ecotype2, geography, ecology, coeffs, linear_model_warning, n_windows_dxy, n_windows_fst)
	row.out$linear_model_warning <- linear_model_warning
	return(row.out)
	
}