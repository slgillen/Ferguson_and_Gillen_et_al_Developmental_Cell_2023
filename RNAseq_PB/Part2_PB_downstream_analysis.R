#this script conducts downstream investigations of the palbociclib RNA-seq data in three neuroblastoma cell lines


# load libraries ----------------------------------------------------------
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(factoextra)
library(scales)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(fgsea)

# collate datasets --------------------------------------------------------

# combine differential expression tables from the different comparisons conducted
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


data_path='PB_RNAseq/output_data/'
filenames <- dir(data_path, pattern = "*.csv")
length(filenames)

allinputs<-lapply(as.list(filenames),function(x) read.delim(paste0(data_path,x),sep=',',stringsAsFactors = FALSE,header=TRUE))
length(allinputs)

allDE_collated<-collate_results(allinputs, filenames)
names(allDE_collated)
dim(allDE_collated)

                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  

