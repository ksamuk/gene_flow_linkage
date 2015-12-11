chrom.to.num <- function(x){
	x <- gsub("group", "", x)
	chrom.rom <- as.character(as.roman(c(1:21)))
	return(match(x, chrom.rom))
}