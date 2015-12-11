####PAML Control file builder
####Builds a control file for a list of phylip files
####requires: paml_file_list.txt, the output of "ls *.phy" in a directory full of phylip files

home.dir<-"~/review/analysis/gma/ev_prep_scripts/paml_analysis"
destination<-"alignments"
setwd(home.dir)

#get the file list and split out the names from the extensions
file.list<-readLines("paml_file_list.txt")
file.list<-sapply(strsplit(file.list,".best.phy"),function(x)x[1])


#seqfile = file.best.phy * sequence data filename
#treefile = file.best.dnd      * tree structure file name
#outfile = file.cml           * main result file name

empty.ctl<-file("empty_control.ctl")
header<-vector(length=3)
for (i in 1:length(file.list)){
  header[1]<-paste("      ","seqfile = ",file.list[i],".best.phy",sep="")
  header[2]<-paste("     ","treefile = ",file.list[i],".best.dnd",sep="")
  header[3]<-paste("      ","outfile = ",file.list[i],".cml",sep="")
  file.name<-paste(file.list[i],".ctl",sep="")
  writeLines(c(header,readLines(empty.ctl)),con=file.path(home.dir,destination,file.name))  
}




