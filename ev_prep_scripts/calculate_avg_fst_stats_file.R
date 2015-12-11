# calculate average fst for every stats file

###############
# LIBRARIES
###############

library("dplyr")
library("parallel")

###############
# FUNCTIONS
###############

calc_average_fst <- function(stats.file.name){
  
  print(paste0("processing ", stats.file.name, "..."))
  
  # read in stats file
  stats.file <- read.table(stats.file.name, header = TRUE, stringsAsFactors = FALSE)
  
  # remove missing, inf and -ve FSTs
  stats.file <- stats.file %>%
    filter(!is.na(Fst)) %>%
    filter(!is.infinite(Fst)) %>%
    filter(Fst > 0)
  
  # calculate average FST (*quotient of FST nums and denoms*, this is the correct way to do this)
  fst.num.sum <- sum(stats.file$FstNum)
  fst.denom.sum <- sum(stats.file$FstDenom)
  fst.avg <- fst.num.sum / fst.denom.sum
  
  # parse the file name
  file.name.stripped <- sapply(strsplit(stats.file.name, split = "/"), function(x)gsub(".txt","",x[length(x)]))
  file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
  pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
  pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
  geography <- file.name.split[3]
  ecology <- file.name.split[4]
  
  # the formatted row of data
  row.out <- data.frame(pop1, ecotype1, pop2, ecotype2, geography, ecology, fst.avg)
  return(row.out)

}

###############
# BODY
###############

stats.file.names <- list.files("stats/snp_all", full.names = TRUE)
#average.fsts <- lapply(stats.file.names, calc_average_fst)
average.fsts <- mclapply(stats.file.names, calc_average_fst, mc.cores = 3, mc.silent = FALSE, mc.preschedule = FALSE) 
average.fsts <- bind_rows(average.fsts)
write.table(average.fsts, file = "meta_data/average_fsts.txt", row.names = FALSE, quote = FALSE)
