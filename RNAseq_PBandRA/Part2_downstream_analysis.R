#this script conducts downstream investigations of the palbociclib and retinoic acid RNA-seq data


# load libraries ----------------------------------------------------------

library(FSA)
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(factoextra)
library(scales)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

# select significantly changed gene set -----------------------------------

collate_results<-function(results_list,results_names){
  for(i in 1:length(results_list)){
    if(i==1){
      all_results<-results_list[[i]][,c('gene','log2FoldChange','padj')]
      names(all_results)[2:3]<-paste0(gsub('_DESeq2output.csv','',results_names[i]),'_',names(all_results)[2:3])
    }else{
      resultsx<-results_list[[i]][,c('gene','log2FoldChange','padj')]
      names(resultsx)[2:3]<-paste0(gsub('_DESeq2output.csv','',results_names[i]),'_',names(resultsx)[2:3])
      all_results<-join(all_results,resultsx,by='gene',type='full')
      rm(resultsx)
    }
  }
  return(all_results)
} 


data_path='PBandRA_RNAseq/outputdata/'
filenames <- dir(data_path, pattern = "*.csv")
length(filenames)

allinputs<-lapply(as.list(filenames),function(x) read.delim(paste0(data_path,x),sep=',',stringsAsFactors = FALSE,header=TRUE))
length(allinputs)

allDE_collated<-collate_results(allinputs, filenames)
names(allDE_collated)
dim(allDE_collated)

# select all genes that have a significant change in at least one condition with a log2FC > 0.5 or < -0.5
alldiffgeneset<-subset(allDE_collated_V2,PBandRAvControl_padj<0.05 | PBandRAvPB_padj<0.05 | PBandRAvRA_padj<0.05 | PBvControl_padj<0.05 | RAvControl_padj<0.05)
nrow(alldiffgeneset) 
alldiffgeneset<-subset(alldiffgeneset, PBandRAvControl_log2FoldChange>0.5 | PBandRAvControl_log2FoldChange<(-0.5) | PBandRAvPB_log2FoldChange>0.5 | PBandRAvPB_log2FoldChange<(-0.5) | PBandRAvRA_log2FoldChange>0.5 | PBandRAvRA_log2FoldChange<(-0.5) | PBvControl_log2FoldChange>0.5 | PBvControl_log2FoldChange<(-0.5) | RAvControl_log2FoldChange>0.5 | RAvControl_log2FoldChange<(-0.5))
nrow(alldiffgeneset) 


# get CPM normalised counts -----------------------------------------------

#raw counts
RNAseq_counts<-read.csv('PBandRA_RNAseq/data/BE2C_PBandRA_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
names(RNAseq_counts)
dim(RNAseq_counts)

#convert to CPM (counts per million)
RNAseq_CPM<-RNAseq_counts
RNAseq_CPM[,2:17]<-apply(RNAseq_CPM[,2:17],2,function(x) (x/sum(x))*1000000)
lapply(RNAseq_CPM[,2:17],function(x) sum(x))

#write output table
write.csv(RNAseq_CPM,file='BE2C_PBandRA_RNAseq_CPM.csv',quote=FALSE,col.names=TRUE,row.names=FALSE)


