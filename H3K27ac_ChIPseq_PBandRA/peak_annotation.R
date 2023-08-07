#!/usr/bin/env Rscript


################ load libraries ###################

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(EnsDb.Hsapiens.v75)




txdb1 <- TxDb.Hsapiens.UCSC.hg19.knownGene

indir<-'' #directory with input files
outdir<-'' #directory for outputs

conditions<-c('RAvDMSO','PBvDMSO','PBRAvDMSO','PBvRA','PBRAvPB','PBRAvRA')
grouptypes<-c('_up05','_ns05','_down05')

########### multi plots ############

plot_ChIPseeker_multiple<-function(txdbX,filenamesX,plotnameX,datadir1,datadir2){
  peakFiles<-lapply(filenamesX, function(filename) readPeakFile(paste0(datadir1,filename,'.bed')))
  names(peakFiles)<-filenamesX
  peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdbX,tssRegion=c(-3000, 3000), verbose=FALSE)
  tiff(paste0(datadir2,'peakAnno_',plotnameX,'.tiff'),height = 10, width = 20, units='cm', compression = "lzw", res = 300)
  print(plotAnnoBar(peakAnnoList))
  dev.off()

  tiff(paste0(datadir2,'distToTSS_',plotnameX,'.tiff'),height = 10, width = 20, units='cm', compression = "lzw", res = 300)
  print(plotDistToTSS(peakAnnoList))
  dev.off()
}

for(cond in conditions){
    filenames<-paste0(cond,grouptypes)
    plot_ChIPseeker_multiple(txdb1,filenames,cond,indir,outdir)
}


########### peak annotations per condition set ############

for(cond in conditions){
    filenames<-paste0(cond,grouptypes)
    for(f in filenames){
        peak <- readPeakFile(paste0(indir,f,'.bed'))
        print(peak)

        txdb1 <- TxDb.Hsapiens.UCSC.hg19.knownGene
        seqlevelsStyle(txdb1) <- "UCSC"
        peakAnno.txdb <- annotatePeak(peak, tssRegion=c(-3000, 3000),TxDb=txdb1, annoDb="org.Hs.eg.db")
        peakAnno.txdb2<-data.frame(peakAnno.txdb,stringsAsFactors=FALSE)
        #peakAnno.txdb2<-apply(peakAnno.txdb2,2,function(x) gsub('\t','',x))
        peakAnno.txdb2$annotation<-gsub('\t','',peakAnno.txdb2$annotation)
        write.table(peakAnno.txdb2,paste0(outdir,'peakAnno_',f,'.txt'),quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')
        rm(peak,txdb1,peakAnno.txdb,peakAnno.txdb2)
    }
}






