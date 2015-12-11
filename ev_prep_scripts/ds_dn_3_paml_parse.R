########parse paml file to data frame
########new approach (uses raw paml files to avoid shell script shenanigans)

library(ape)
library(biomaRt)
library(dplyr)
library(magrittr)

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/review/analysis/gma/ev_prep_scripts/paml_analysis"

#the output directory (window evs)
out.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/window"

#contains pre-processed paml output (the dn and ds tree lines+first file of paml file)
paml.output.dir<-file.path(home.dir,"alignments_all")
setwd(paml.output.dir)

#dat file list
file.list<-list.files()

#prep output data frame
gene.id<-vector(mode="character",length=length(file.list))
ds<-vector(mode="numeric",length=length(file.list))
dn<-vector(mode="numeric",length=length(file.list))
num.sites<-vector(mode="numeric",length=length(file.list))
gacu.dnds<-data.frame(gene.id,ds,dn,num.sites,stringsAsFactors=FALSE)

#loop through files and extract ds info
for (i in 1:length(file.list)){
  
  #read in file
  file.con<-file(file.list[i])
  file.lines<-readLines(file.con)
  closeAllConnections()
  
  #find the "ds tree" line
  ds.tree.line<-grep("dS tree",file.lines)
  
  #find the "After deleting gaps. " line (for filtering bad alignments)
  gaps.line<-grep("After deleting gaps. ",file.lines)
  num.sites<-as.numeric(unlist(strsplit(file.lines[gaps.line],split=" "))[4])
  
  #edge case where there are no alignment gaps
  if(length(num.sites)==0){
    num.sites<-as.numeric(unlist(strsplit(file.lines[1],split=" ")))
    num.sites<-num.sites[length(num.sites)]
  }
  gacu.dnds$num.sites[i]<-num.sites
    
  #find gene name from the file name
  gene.name<-sapply(strsplit(file.list[i],split=".cml"),function(x)x[1])
  if(gene.name==""){
    print(file.list[i]," has no gene name?!")
  }
  gacu.dnds$gene.id[i]<-gene.name
  
  #if no ds tree or target gene, skip file
  #uncomment for commentary on the quality of your data files
  if(length(ds.tree.line)==0){
    #print(paste(file.list[i],"is missing dS tree."))
    gacu.dnds$ds[i]<-NA
    gacu.dnds$dn[i]<-NA
  }else if(length(grep(paste(gene.name,":",sep=""),file.lines))==0){
    #print(paste(file.list[i],"has a dS tree, but is missing the target gene."))
    gacu.dnds$ds[i]<-NA
    gacu.dnds$dn[i]<-NA
  }else{
    
    #make the tree(s)
    ds.tree<-read.tree(text=file.lines[ds.tree.line+1])
    dn.tree<-read.tree(text=file.lines[ds.tree.line+3])
    
    #if there is no dn or ds value, assign NA, otherwise grab value from tree
    if(is.null(ds.tree$edge.length[which.edge(ds.tree,gene.name)])){
      gacu.dnds$ds[i]<-NA
    }else{
      gacu.dnds$ds[i]<-ds.tree$edge.length[which.edge(ds.tree,gene.name)]
    }
    
    if(is.null(dn.tree$edge.length[which.edge(dn.tree,gene.name)])){
      gacu.dnds$dn[i]<-NA
    }else{
      gacu.dnds$dn[i]<-dn.tree$edge.length[which.edge(dn.tree,gene.name)]
    }
    
    #progress bar
    if (i%%1000==0){
      cat(i,"...",sep="")
    }
  }
}

#remove lines containing NAs
gacu.dnds<-gacu.dnds[complete.cases(gacu.dnds),]

##match gene.ids to genomic coordinates

#initialize gacu ensembl
ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")

#the attributes of interest
attributes.feat<-c("ensembl_gene_id",
                   "ensembl_peptide_id",
                   "start_position",
                   "end_position",
                   "chromosome_name")

#query ensembl for coordinates
coords<-getBM(attributes=attributes.feat,values=gacu.dnds$gene.id,filters=c("ensembl_peptide_id"),mart = ensembl)

#for easier viewing
coords<-arrange(coords,ensembl_peptide_id)

#match in ds/dn values (could be a left_join, whateves)
ds.out<-gacu.dnds$ds[match(coords$ensembl_peptide_id,gacu.dnds$gene.id)]
dn.out<-gacu.dnds$dn[match(coords$ensembl_peptide_id,gacu.dnds$gene.id)]

#build prelim output file
gacu.out<-data.frame(gene.id=coords$ensembl_gene_id,
                     peptide.id=coords$ensembl_peptide_id,
                     lg=coords$chromosome_name,
                     pos1=coords$start_position,
                     pos2=coords$end_position,
                     ds=ds.out,
                     dn=dn.out,
                     sites=num.sites)

#if there are multiple peptides from a single gene (~30% of data), take the mean of the ds values
out.means<-gacu.out%>%
  group_by(gene.id)%>%
  summarise(mean(ds),mean(dn))
out.means<-data.frame(out.means)
names(out.means)<-c("gene.id","ds.mean","dn.mean")

#join in means with rest of data
gacu.out.2<-left_join(gacu.out,out.means)

#remove scaffolds
gacu.out.2<-gacu.out[grep("group",gacu.out$lg),]
#convert lg to numeric
gacu.out.2$lg<-as.numeric(as.roman(sapply(strsplit(as.character(gacu.out.2$lg),"group"),function(x)x[2])))
#arrange
gacu.out.2<-arrange(gacu.out.2,lg)


#FILTERING: no dn/dsvalues above 2 (from literature), no genes with <=100 sites
gacu.out.2<-gacu.out.2[gacu.out.2$sites>=100,]
gacu.out.2<-gacu.out.2[gacu.out.2$dn<=2,]
gacu.out.2<-gacu.out.2[gacu.out.2$ds<=2,]

#prep output files (strip duplicates)
gacu.out.ds<-gacu.out.2[,3:6]
gacu.out.ds<-unique(gacu.out.ds)

gacu.out.dn<-gacu.out.2[,c(3:5,7)]
gacu.out.dn<-unique(gacu.out.dn)


#write to file
setwd(out.dir)
write.table(gacu.out.ds,file="ds.txt",row.names=FALSE)
write.table(gacu.out.dn,file="dn.txt",row.names=FALSE)



