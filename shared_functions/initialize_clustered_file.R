########################################
# Initialize cluster output file
# (used by all cluster permutation scripts)
########################################

initialize_clustering_output <- function(){
	
## the clustering file
cluster.df <- read.table(file = "analysis_ready/snp_clustering_metrics.txt", header = TRUE, stringsAsFactors = FALSE)
cluster.df <- cluster.df %>%
	mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp) %>%
	mutate(nnd.diff.sd = nnd.diff/nnd.sd.null)

## add in region data (for looser geography)
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short pop codes
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

# associate pops in coeff dat with regions in region dat
region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop1","reg1")
cluster.df$reg1 <- region.sub$reg1[match(cluster.df$pop1, region.sub$pop1)]

region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop2","reg2")
cluster.df$reg2 <- region.sub$reg2[match(cluster.df$pop2, region.sub$pop2)]

# make new geographic categories
cluster.df$geography2 <- ifelse(cluster.df$reg1==cluster.df$reg2, "para", "allo")

#make new groups 
cluster.df$group2 <- paste0(cluster.df$geography2,"_",cluster.df$ecology)

group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
group.rename <- c("Allopatry\nDivergent", "Allopatry\nParallel", "Gene Flow\nDivergent", "Gene Flow\nParallel")
cluster.df$group2.new <-group.rename[match(cluster.df$group2, group.old.names)]
cluster.df$group.new <- group.rename[match(cluster.df$group, group.old.names)]
return(cluster.df)
}