## add in region data (for looser geography)
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short pop codes
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

parse_stats_file_name <- function (stats.file.name){
  
  file.name.stripped <- sapply(strsplit(stats.file.name, split = "/"), function(x)gsub(".txt","",x[length(x)]))
  file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
  pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
  pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
  geography <- file.name.split[3]
  ecology <- file.name.split[4]
  group <- paste0(geography, "_", ecology)
  comparison <- file.name.stripped
  row.out <- data.frame(pop1, ecotype1, pop2, ecotype2, geography, ecology, comparison,group)

  
  # associate pops in coeff dat with regions in region dat
  region.sub <- region.dat[,c(1,4)]
  names(region.sub) <- c("pop1","reg1")
  row.out$reg1 <- region.sub$reg1[match(row.out$pop1, region.sub$pop1)]
  
  region.sub <- region.dat[,c(1,4)]
  names(region.sub) <- c("pop2","reg2")
  row.out$reg2 <- region.sub$reg2[match(row.out$pop2, region.sub$pop2)]
  
  # make new geographic categories
  row.out$geography2 <- ifelse(row.out$reg1 == row.out$reg2, "para", "allo")
  
  #make new groups :o
  row.out$group2 <- paste0(row.out$geography2, "_", row.out$ecology)
  
  return(row.out)

}