####gene metadata masterfile
####binds ensembl metadata with a list of genomic coordinates (from FST/stats files)
####outputs ev file for each variable of interest (must be one per file to work with ev binding script)
####INPUTS: ensembl_exon_positions.txt, ensembl_gene_metadata, all_stats_all_pops_date.csv
####OUTPUTS: gene_metadata.txt

#libraries
library("dplyr")
library("fastmatch")

###input files
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis"
output.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/window"
setwd(home.dir)

#grab the positions from an all_stats file
stats.positions<-read.csv(file="all_stats_all_pops_feb9-2014.csv",header=TRUE)
stats.positions<-stats.positions[,1:2]

setwd(file.path(home.dir,"ev prep scripts"))

exon.positions<-read.table(file="ensembl_exon_positions.txt",header=TRUE)
gene.positions<-read.table(file="ensembl_gene_metadata.txt",header=TRUE)

###trim input files of scaffolds
gene.positions<-gene.positions[grep("group",gene.positions$chromosome_name),]
exon.positions<-exon.positions[exon.positions$ensembl_gene_id%in%gene.positions$ensembl_gene_id,]

#chrom to numeric
gene.positions$chromosome_name<-as.numeric(as.roman(sapply(strsplit(as.character(gene.positions$chromosome_name),split="group"),function(x)x[2])))
#naming to match ev script
names(gene.positions)[2:4]<-c("lg","pos1","pos2")

#utrs
utr5.positions<-data.frame(lg=gene.positions$lg)
utr5.positions$pos1<-sapply(strsplit(as.character(gene.positions$utr5),split=":"),function(x)x[1])
utr5.positions$pos2<-sapply(strsplit(as.character(gene.positions$utr5),split=":"),function(x)x[2])
utr5.positions$utr5<-TRUE

#utr3
utr3.positions<-data.frame(lg=gene.positions$lg)
utr3.positions$pos1<-sapply(strsplit(as.character(gene.positions$utr3),split=":"),function(x)x[1])
utr3.positions$pos2<-sapply(strsplit(as.character(gene.positions$utr3),split=":"),function(x)x[2])
utr3.positions$utr3<-TRUE

#get lgs for the exon positions
#fastmatch! its fast.
detach("package:fastmatch", unload=TRUE)
exon.positions$lg<-gene.positions$lg[match(exon.positions$ensembl_gene_id,gene.positions$ensembl_gene_id)]

#split exons into pos ranges
exon.positions$pos1<-sapply(strsplit(as.character(exon.positions$exon),split=":"),function(x)x[1])
exon.positions$pos2<-sapply(strsplit(as.character(exon.positions$exon),split=":"),function(x)x[2])
exon.positions$exon<-TRUE

####BUILD OUTPUT FILES

#ensembl_gene_id (gene y/n)
gene.positions.gene_id<-gene.positions%>%
  select(lg,pos1,pos2,ensembl_gene_id)
#exon positons
exon.positions<-exon.positions%>%
  select(lg,pos1,pos2,exon)
#utr positions previously built

#write output files to disk
setwd(output.dir)
write.table(gene.positions.gene_id,file="gene_id.txt",row.names=FALSE)
write.table(exon.positions,file="exon.txt",row.names=FALSE)
write.table(utr5.positions,file="utr5.txt",row.names=FALSE)
write.table(utr3.positions,file="utr3.txt",row.names=FALSE)
