#!/usr/bin/env Rscript

files <- list.files(pattern = "\\.gff$")
length(files)

for(f in files){
  xx<-read.delim(f,header=FALSE,stringsAsFactors = FALSE,sep='\t')
  xx$size<-xx$V5-xx$V4
  write.table(xx[,c(1,4,5,2,10,7)],file=gsub('.gff','.bed',f),sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
  rm(xx)
}
