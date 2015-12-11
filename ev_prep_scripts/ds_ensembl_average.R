#### Alternate dS approach: the mean of all ensembl dS estimates


#installing biomaRt from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(lib="~/bin/R/libraries")
biocLite(c("biomaRt","rtracklayer"),lib="~/bin/R/libraries")
biocLite(c("rtracklayer"))
#libraries
library(biomaRt)
library(dplyr)
library(data.table)

#set the mart
ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")
filters <- listFilters(ensembl)
attributes <-  listAttributes(ensembl, what=c("name","description","page"))

# the sets of attributes to be obtained from ensembl
attributes.feat <- c("ensembl_gene_id",
                     "start_position",
                     "end_position",
                     "chromosome_name")

attributes.struct <- c("ensembl_gene_id",
                       "start_position",
                       "end_position",
                       "genomic_coding_start",
                       "genomic_coding_end")

# all the attributes with "_ds" in their name
ds.homologs <- attributes %>%
  filter(page == "homologs") %>%
  filter(grepl("_ds",name)) %>%
  select(name) %>%
  unlist(use.names = FALSE)

# query ensembl
data.feat <- getBM(attributes=attributes.feat,values="*",mart = ensembl)
data.struct <- getBM(attributes=attributes.struct,values="*",mart = ensembl)
data.homologs <- getBM(attributes=ds.homologs[1:6],values="*",mart = ensembl)


gene.ids<- "ENSGACG00000018311"

gene.ids <- data.feat$ensembl_gene_id[10000]

data.homologs <- getBM(attributes=ds.homologs[46], filter="ensembl_gene_id", values="gene.ids" , mart = ensembl)

calculate_ds <- function(gene.id) {
  
}

#unique gene ids from above two queries, fed to seq query below
master.gene.ids<-unique(c(master.list.struct$ensembl_gene_id,master.list.feat$ensembl_gene_id))
#master.gene.ids<-master.gene.ids[1:10]
#master.list.seq<-getBM(attributes=master.attributes.seq,filters=c("ensembl_gene_id"),values=master.gene.ids,mart = ensembl)

#master.list.seq<-getBM(attributes=master.attributes.seq,filters=c("ensembl_gene_id"),values="ENSGACG00000000002",mart = ensembl)

#combine duplicate rows in feature dataset
#(geneIDs duplicated if gene has >1 GO term)

DT <- data.table(master.list.feat, key="ensembl_gene_id")
DT2<-DT[, list(chromosome_name=chromosome_name,
               start_position=start_position, 
               end_position=end_position, 
               gene_biotype=gene_biotype,
               go_ids=list(go_id)),
        by=ensembl_gene_id]

go.df<-data.frame(ensembl_gene_id=DT$ensembl_gene_id,go_id=DT$go_id)

#collapse duplicated rows
master.list.feat.collapsed<-data.frame(DT2)
master.list.feat.collapsed<-master.list.feat.collapsed[!duplicated(master.list.feat.collapsed$ensembl_gene_id),]
