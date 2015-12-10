#### process raw clustering script output

################################################################################
# Libraries and initalizing variables
################################################################################
library("dplyr")

################################################################################
# Input files
################################################################################
cluster.folder <- "analysis_ready/clustering_fst_new"
cluster.files <- list.files(cluster.folder, full.names = TRUE)

# grab the header from only the first file
header.row <- names (read.table(file = cluster.files[1], header = TRUE, stringsAsFactors = FALSE))

################################################################################
# Process raw output
################################################################################

# read in files w/o headers
read_files_chopped <- function(cluster.file){
  return(read.table(file = cluster.file, stringsAsFactors = FALSE, skip = 1))
}

file.chopped <- lapply(cluster.files, read_files_chopped)
cluster.df <- bind_rows(file.chopped)
names(cluster.df) <-  header.row

# meta data for relaxed groups
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short population codes from file names
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

# associate pops in cluster df with regions in region dat
region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop1","reg1")
cluster.df$reg1 <- region.sub$reg1[match(cluster.df$pop1, region.sub$pop1)]

region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop2","reg2")
cluster.df$reg2 <- region.sub$reg2[match(cluster.df$pop2, region.sub$pop2)]

# make relaxed geographic categories
cluster.df$geography2 <- ifelse(cluster.df$reg1==cluster.df$reg2, "para", "allo")

# make relaxed groups
cluster.df$group2 <- paste0(cluster.df$geography2,"_",cluster.df$ecology)

# arrange columns
cluster.df <- cluster.df %>%
  mutate(group = paste0(geography, "_", ecology)) %>%
  select(pop1, ecotype1, reg1, pop2, ecotype2, reg2, geography, geography2, ecology, group, group2, everything())

# kill duplicates
cluster.df <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(comparison, lg) %>%
	filter(row_number(comparison) == 1)


################################################################################
# write to file
################################################################################

write.table(cluster.df, file = "analysis_ready/snp_clustering_metrics.txt", quote = FALSE, row.names = FALSE)
