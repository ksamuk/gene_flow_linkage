plot_averaged_regression_lines <- function(coeff.dat, group_order, pal, group_variable, stat = "fst", ylim = c(0,1), xlim = c(0,1), ylab = ylab, adjust = 0.0005, label = ""){
	
	recomb_rate_type <- paste0("recomb_rate_", stat)
	intercept_type <- paste0("intercept_", stat)
	
	coeff.dat <- coeff.dat %>% 
		mutate_(recomb_rate = interp(recomb_rate_type)) %>%
		mutate_(intercept = interp(intercept_type))
	
	# set up the plot
	label.scaling <-  1.5
	line.weight <- 12
	par(mgp = c(2.5, 1, 0), mar = c(5, 4, 2, 1), mex = 1, bty = "l")
	plot(1, ylim = ylim, xlim = xlim, 
			 xlab = "Recombination Rate (cM/Mb)", ylab = ylab, 
			 cex.lab = label.scaling, cex.axis = label.scaling, yaxs= "i", xaxs = "i")
	text(49, (ylim - ylim/80), label, cex = label.scaling)
	
	# plot each line individually
	for (i in 1:length(group_order)){
		
		# subset data by group
		filter_criteria <- interp(~ which_column == group_order[i], which_column = as.name(group_variable))
		coeff.dat.group <- coeff.dat %>%
			filter_(filter_criteria)
		
		# average intcepts and slopes
		intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
		slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
		
		#define an equation function
		#0.0005 is for style only (keeps line from overplotting axis)
		eq <- function(x){
			ylog <- intercept + slope*x 
			return(1 - 1/(1 + exp(ylog)) + adjust)
		} 
		
		col.line <- pal[i]
		curve(eq, from = 0, to = 50, ylim = ylim, add=TRUE, lwd = line.weight, col = col.line)
		
	}
	
}