#builds metadata/ev file from biomart
#using biomaRt library

#working directory setup
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev prep scripts/ds_processing"
#setwd(home.dir)
#rm(list=ls())

#installing biomaRt from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(lib="~/bin/R/libraries")
biocLite(c("biomaRt","rtracklayer"),lib="~/bin/R/libraries")
biocLite(c("rtracklayer"))
#libraries
library(rtracklayer)
library(biomaRt)
library(dplyr)
library(data.table)

#set the mart
ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")
#filters <- listFilters(ensembl)
attributes <-  listAttributes(ensembl, what=c("name","description","page"))

#query emsembl for master list of protein coding genes

#attributes to pull from ensembl
#feat(ure), struct(ure), seq(quences) are "pages" of ensembl attributes 
#can only access attrs from one page at a time, so three separate queries are needed
master.attributes.feat<-c("ensembl_gene_id",
                     "start_position",
                     "end_position",
                     "chromosome_name",
                     "gene_biotype",
                     "go_id",
                     "name_1006",
                     "percentage_gc_content")

master.attributes.struct<-c("ensembl_gene_id",
                            "start_position",
                            "end_position",
                            "genomic_coding_start",
                            "genomic_coding_end",
                            "5_utr_start",
                            "5_utr_end",
                            "3_utr_start",
                            "3_utr_end")

master.list.feat<-getBM(attributes=master.attributes.feat,values="*",mart = ensembl)
master.list.struct<-getBM(attributes=master.attributes.struct,values="*",mart = ensembl)

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

#combine duplicate rows in structure dataset
#(geneIDs duplicated if gene has >1 exon)

#R doesn't abide column names starting with numerals
names(master.list.struct)[6:9]<-c("utr5_start","utr5_end","utr3_start","utr3_end")

DT <- data.table(master.list.struct, key="ensembl_gene_id")
DT2<-DT[, list("genomic_coding_start"=paste(genomic_coding_start,genomic_coding_end,sep=":"),
               "utr5"=paste(utr5_start,utr5_end,sep=":"),
               "utr3"=paste(utr3_start ,utr3_end ,sep=":")),by=ensembl_gene_id]

DT3<-DT2[, list(exons=list(genomic_coding_start),
               utr5=utr5,
               utr3=utr3),
        by=ensembl_gene_id]

#grab exons (must be written separately to avoid dups and writing lists to file)
exons.df<-data.frame(ensembl_gene_id=DT2$ensembl_gene_id,exon=DT2$genomic_coding_start))

#excise duplicates
master.list.struct.collapsed<-data.frame(DT3)
master.list.struct.collapsed<-master.list.struct.collapsed[!duplicated(master.list.struct.collapsed$ensembl_gene_id),]

#merged metadata file
tmp<-merge(master.list.feat.collapsed,master.list.struct.collapsed,by="ensembl_gene_id")

#cut out go ids and exons (now written separately)
tmp.trim<-tmp[,-c(6,7)]

#write metadata to file
write.table(tmp.trim,file="ensembl_gene_metadata.txt")

#write exons to file
write.table(exons.df,file="ensembl_exon_positions.txt")

#write individual go terms to file (grabbed before feat DT)
write.table(go.df,file="ensembl_go_terms.txt")

#extract master go term definitions and write to file
go.out<-data.frame(go_id=master.list.feat$go_id,go_cat=master.list.feat$name_1006)
go.out<-go.out[unique(go.out$go_id),]
write.table(go.out,file="all_go_ids_categories_ensembl.txt")


