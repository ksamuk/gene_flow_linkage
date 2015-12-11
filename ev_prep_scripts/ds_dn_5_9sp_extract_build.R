#pulling ds values from ninespine data
#goal is to compare ks values from baocheng to our paml results

library(dplyr)
library(ggplot2)
library(biomaRt)

working.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/raw"
setwd(working.dir)

kaks<-read.table("one2one_kaks_split.txt",header=TRUE)

# define biomart object
ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")
master.attributes.feat<-c("ensembl_gene_id",
                          "ensembl_peptide_id",
                          "ensembl_transcript_id",
                          "start_position",
                          "end_position",
                          "chromosome_name")

#grab attributes for peptide IDs from kaks file
master.list.feat<-getBM(attributes=master.attributes.feat,
                        filters = "ensembl_peptide_id", 
                        values = unique(kaks$X3spgene),
                        mart = ensembl)

#match in ka/ks values
master.list.feat$ks<-kaks[match(master.list.feat$ensembl_peptide_id,kaks$X3spgene),]$Ks
master.list.feat$ka<-kaks[match(master.list.feat$ensembl_peptide_id,kaks$X3spgene),]$Ka

#only genes from lgs
master.list.out<-master.list.feat[grep("group",master.list.feat$chromosome_name),]

#lg to numeric
master.list.out$chromosome_name<-as.numeric(as.roman(sapply(strsplit(master.list.out$chromosome_name,split="group"),function(x)x[2])))

#build output files
ks.out<-data.frame(lg=master.list.out$chromosome_name,
                   pos1=master.list.out$start_position,
                   pos2=master.list.out$end_position,
                   ks=master.list.out$ks)
ks.out<-arrange(ks.out,lg,pos1)

ka.out<-data.frame(lg=master.list.out$chromosome_name,
                   pos1=master.list.out$start_position,
                   pos2=master.list.out$end_position,
                   ka=master.list.out$ka)
ka.out<-arrange(ka.out,lg,pos1)

#write to file

out.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/window"
setwd(out.dir)

write.table(ks.out,file="ks.txt",row.names=FALSE)
write.table(ka.out,file="ka.txt",row.names=FALSE)
