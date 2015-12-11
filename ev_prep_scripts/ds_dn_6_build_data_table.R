####create data table from (pre-processed) codeml output
####outputs gene ids and ds values for gacu ONLY
####KS mar-3-2015

library(ape)
library(biomaRt)
library(dplyr)

# as if this wasn't going to be here
chrom.to.num <- function(x){
  x <- gsub("group", "", x)
  chrom.rom <- as.character(as.roman(c(1:21)))
  return(match(x, chrom.rom))
}

#home dir set up
#home.dir<-"E:/Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"

#contained pre-processed paml output (the dn and ds trees)
paml.output.dir <- file.path("ev_prep_scripts/paml/output")

#dat file list
file.list <- list.files(paml.output.dir, full.names = TRUE)

#set up empty output dataframe
gene.id <- vector(mode = "character", length = length(file.list))
ds <- vector(mode = "numeric", length = length(file.list))
gacu.ds <- data.frame(gene.id, ds, stringsAsFactors = FALSE)

for (i in 1:length(file.list)){
  #read in file
  file.con <- file(file.list[i])
  file.lines <- readLines(file.con)
  close(file.con)
  if (!is.na(file.lines[2])){  
    # make the tree
    ds.tree <- read.tree(text = file.lines[3])
    
    # distance matrix for branches
    dist.mat <- cophenetic.phylo(ds.tree)
    
    # the gacu row in the distance matrix
    gacu.row <- row.names(dist.mat) %>% grep("ENSGACP",.) %>% dist.mat[.,]
    
    # only retain estimates for pformosa, oniloticus, xmaculatus
    gacu.row <- gacu.row %>% names %>% grep("pformosa|oniloticus|xmaculatus",.) %>% gacu.row[.]
    
    # remove implausible branch lengths (likely bad alignments)
    # the value of "3" was chosen because it was the 99th percentile of the distribution of estimated ds values
    gacu.row <- gacu.row[gacu.row < 3]
    
    # the ds estimate 
    gacu.ds$ds[i] <- mean(gacu.row, na.rm = TRUE)
    
    #get gene name and find it in the ds tree
    gene.name <- strsplit(file.list[i], split=".cml.txt") %>% gsub(paml.output.dir,"",.) %>% gsub("/","",.)
    
    #pull out ds value for gacu
    gacu.ds$gene.id[i] <- gene.name
    
    #tracker
    if (i %% 1000 == 0){
      print(i)
    }
    
  }
  if (is.na(file.lines[1])){
    gacu.ds$gene.id[i] <- NA
    gacu.ds$ds[i] <- NA
  }
  
}

# add in locations

# the physical positons for all gacu genes
#gene.dat <- read.table(file = "evs/additional/gene_id.txt", stringsAsFactors = FALSE, header = TRUE)

# convert prot ids to gene ids (???)

ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")
master.attributes.struct<-c("ensembl_gene_id",
                            "ensembl_peptide_id",
                            "start_position",
                            "end_position", 
                            "chromosome_name")

master.list.struct <- getBM(attributes = master.attributes.struct, filters = "ensembl_peptide_id", values = gacu.ds$gene.id, mart = ensembl)

names(gacu.ds)[1] <- "ensembl_peptide_id" 

gacu.ds <- left_join(master.list.struct, gacu.ds[,1:2], by = "ensembl_peptide_id")

gacu.ds <- gacu.ds %>%
  filter(!is.nan(ds)) %>%
  filter(grepl("group",chromosome_name)) %>%
  group_by(ensembl_gene_id) %>%
  summarise(lg = chromosome_name, pos1 = mean(start_position), pos2 = mean(end_position), ds = mean(ds)) %>%
  select(lg, pos1, pos2, ds)

gacu.ds$lg <- chrom.to.num(gacu.ds$lg)

gacu.ds <- gacu.ds %>%
  arrange(lg)

write.table(gacu.ds, file="evs/ds.txt")

