###building fasta files for homologous genes
###will be used to manually calculate ds using PAML
###3p + sp +10 other teleosts
###pulls from bioMart

library(dplyr)
library(ggplot2)
library(biomaRt)

working.dir<-"/Users/Kieran/Documents/Science/Projects/Ph.D./Genome\ Meta\ Analysis/ev\ prep\ scripts"
setwd(working.dir)

#ensembl<-read.table("ensembl_gene_metadata.txt",header=TRUE)

#1:1 homologous gene pairs for 3ps and 9sp 
#from Baocheng Guo via email
kaks<-read.table("one2one_kaks_split.txt",header=TRUE)

# define biomart object
mart <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")
attributes <-  listAttributes(mart, what=c("name","description","page"))

filters <-  listFilters(mart, what=c("name","description"))

##note: ENSGACP and ENSGACG do not map 1:1
#convert to geneids?
#actually will leaves as peptide ids and convert to transcript ids at the end

gene.ids <- getBM(attributes = c("ensembl_peptide_id","ensembl_transcript_id"),
                  filters = "ensembl_peptide_id", values = unique(kaks$X3spgene),
                  mart = mart)

#query biomart for homologous sequences
# takifugu rubripes (fugu), tetraodon nigroviridis (puffer), oniloticus (nile tilapia), xmaculatus (platyfish) pformosa (amazon molly)

#list of target species
homology.peptide.id<-c("pformosa_homolog_ensembl_peptide",
                  "oniloticus_homolog_ensembl_peptide",
                  "xmaculatus_homolog_ensembl_peptide",
                  "tnigroviridis_homolog_ensembl_peptide",
                  "trubripes_homolog_ensembl_peptide")

homology.ortho.type<-c("pformosa_homolog_orthology_type", 
                  "oniloticus_homolog_orthology_type", 
                  "xmaculatus_homolog_orthology_type", 
                  "tnigroviridis_homolog_orthology_type", 
                  "trubripes_homolog_orthology_type")

#find all the one2one homologs based on the above list
#this might seem like an odd approach, but ensembl does not let you selectively download one2one orthologs
#so, a loop 

homologs<-list()
for (i in 1:length(homology.gene.id)){
  homologs.tmp <- getBM(attributes = c("ensembl_transcript_id",homology.peptide.id[i],homology.ortho.type[i]),
                    filters = "ensembl_transcript_id", values = unique(gene.ids$ensembl_transcript_id),
                    mart = mart)

  homologs[[i]]<-subset(homologs.tmp[,1:2],homologs.tmp[3]=="ortholog_one2one")
  
}

#merge into a single dataframe

merge.all <- function(...){
  merge(..., all=TRUE,by="ensembl_transcript_id")
}

out <- Reduce(merge.all, homologs)

#query each ensembl genome using the list of orthologs
#download sequences into a list
#included gacu for checking purposes

data.sets <- c("pformosa_gene_ensembl",
           "oniloticus_gene_ensembl",
           "xmaculatus_gene_ensembl",
           "tnigroviridis_gene_ensembl",
           "trubripes_gene_ensembl")

gene_exon_intron
gene_exon

for (i in 1:length(data.sets)){
  mart <- useMart("ensembl",dataset=data.sets[1])
  sequences.tmp <- getBM(attributes = c("ensembl_peptide_id","cdna"),
                        filters = "ensembl_peptide_id", values = homologs[[1]][,2],
                        mart = mart)
  homologs[[1]][,3]<-sequences.tmp$ensembl_peptide_id
}


