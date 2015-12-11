### extract a linkage group from a stats file and format it for nnd calculation


extract_lg <- function (file){
  name.split <- strsplit(file,split = "/") 
  name.split <- name.split[[1]][3]
  name.split <- strsplit(name.split,split = "[.]") %>% unlist
  comparison <- paste0(name.split[1],".",name.split[2])
  group <- paste0(name.split[3],"_",name.split[4])
  dat.tmp <- read.table(file = file, header = TRUE, stringsAsFactors = FALSE)
  
  dat.tmp <- dat.tmp %>%
    filter(CHROM == "groupIV") %>%
    filter(!is.infinite(Fst), !is.na(Fst), Fst > 0) %>%
    select(CHROM, POS, Fst)
  
  names(dat.tmp)[1:3] <- c("lg","pos","fst")
  dat.tmp$comparison <- comparison
  dat.tmp$group <- group
  dat.tmp$lg <- 4
  dat.tmp$fst.outlier <- is.outlier(dat.tmp$fst)
  dat.tmp <- dat.tmp %>%
    filter(fst.outlier == TRUE)
  
  #add map distances (generates some warnings, but works as intended)
  dat.tmp <- add_map_distance(dat.tmp)
  return(dat.tmp)
}