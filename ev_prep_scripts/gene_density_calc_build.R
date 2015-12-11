####Gene density calulator
####Outputs two evs: "gene_density" and "gene_coverage"
####gene_density = number of genes in the window
####    _coverage = proportion of window that is "gene"

#washy washy
rm(list=ls())

library(rtracklayer)
library(biomaRt)
library(dplyr)
library(data.table)

########BIOMART: GENE LOCATIONS

#set the biomart
gacu.ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")

#gene attributes from the "feature" page
attributes.feat<-c("start_position",
                    "end_position",
                    "chromosome_name")

gene.list.raw<-getBM(attributes=attributes.feat,values="*",mart = gacu.ensembl)

#remove scaffolds
gene.list<-gene.list.raw[grep("group",gene.list.raw$chromosome_name),]

#chromo to numeric
gene.list$chromosome_name<-as.numeric(as.roman(gsub("group","",gene.list$chromosome_name)))
gene.list<-arrange(gene.list,chromosome_name,start_position)

#standard formatting
gene.list<-data.frame(lg=gene.list$chromosome_name,pos1=gene.list$start_position,pos2=gene.list$end_position)

########END BIOMART: GENE LOCATIONS

########SLIDING WINDOWS OF GENE COUNT

window.size<-75000

#kind of a grody implementation

gene.counts.chr<-list()
for (i in 1:max(gene.list$lg)){
  gene.list.chr<-subset(gene.list,gene.list$lg==i)
  
  #convert nuc position to window index
  window.nums<-as.integer(gene.list.chr$pos1/window.size)+1
  
  #tally up genes in windows
  window.names<-as.numeric(names(table(window.nums)))
  window.counts<-as.numeric(table(window.nums))
  
  #reformatting 
  gene.counts.chr[[i]]<-data.frame(win.num=window.names,gene.count=window.counts)
  
  #calc window boundaries
  window.pos1<-(gene.counts.chr[[i]]$win.num-1)*window.size
  window.pos1[window.pos1==0]<-1
  window.pos2<-((gene.counts.chr[[i]]$win.num)*window.size)-1
  
  #the output for this chromo
  lg<-rep(i,length(window.pos2))
  gene.counts.chr[[i]]<-data.frame(lg=lg,pos1=window.pos1,pos2=window.pos2,gene.count=window.counts)
}

#bind
gene.counts.out<-do.call("rbind",gene.counts.chr)

#write to file
write.table(gene.counts.out,file=file.path(getwd(),"evs","gene_count.txt"),row.names=FALSE)

########END SLIDING WINDOWS OF GENE COUNT

########SLIDING WINDOWS OF GENE DENSITY

#build fake "density" measure
gene.list.out<-gene.list
gene.list.out$gene.density<-1

prop.overlap<-list()
for (i in 1:max(gene.list.out$lg)){

  gene.dens.chr<-subset(gene.list.out,gene.list.out$lg==i)
  windows.chr<-subset(gene.counts.out,gene.counts.out$lg==i)

  ##from "match_evs_to_outlier_file"
  
  ev.range<-IRanges(start=gene.dens.chr$pos1,end=gene.dens.chr$pos2)
  ev.range<-reduce(ev.range)
  stat.range<-IRanges(start=windows.chr$pos1,end=windows.chr$pos2)
  
  #find ovelaps amd build an "overlap df"
  overlap<-findOverlaps(stat.range,ev.range,select="all")
  overlap.df<-data.frame(lg=windows.chr[queryHits(overlap),]$lg,
                         pos1=windows.chr[queryHits(overlap),]$pos1,
                         pos2=windows.chr[queryHits(overlap),]$pos2,
                         gene.start=start(ev.range[subjectHits(overlap)]),
                         gene.end=end(ev.range[subjectHits(overlap)]))
  
  #truncate overlaps that start before window
  overlap.df$gene.start[overlap.df$gene.start<overlap.df$pos1]<-overlap.df$pos1[overlap.df$gene.start<overlap.df$pos1]
  
  #truncate overlaps that extend past end of window
  overlap.df$gene.end[overlap.df$gene.end>overlap.df$pos2]<-overlap.df$pos2[overlap.df$gene.end>overlap.df$pos2]
  
  #recalc widths
  overlap.df$width<-(overlap.df$gene.end-overlap.df$gene.start)+1
  
  #calculate the proportion of total overlap
  prop.overlap[[i]]<-overlap.df%>%
    group_by(pos1)%>%
    summarise(gene.density=sum(width)/window.size)
  
  prop.overlap[[i]]<-data.frame(prop.overlap[[i]])
  prop.overlap[[i]]<-left_join(prop.overlap[[i]],windows.chr)
  prop.overlap[[i]]<-select(prop.overlap[[i]],lg,pos1,pos2,gene.density)

}

gene.dens.out<-do.call("rbind",prop.overlap)

head(gene.dens.out)
tail(gene.dens.out)

#write to file
write.table(gene.dens.out,file=file.path(getwd(),"evs","gene_density.txt"),row.names=FALSE)

########END SLIDING WINDOWS OF GENE DENSITY
