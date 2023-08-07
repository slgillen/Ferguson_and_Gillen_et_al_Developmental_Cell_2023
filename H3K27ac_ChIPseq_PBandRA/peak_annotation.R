#!/usr/bin/env Rscript


# load libraries ------------------------------------------

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)


# plot peak annotation information for H3K27ac response groups ------------------------------------------
indir<-'' #input directory
outdir<-'' #output directory

plot_ChIPseeker_multiple<-function(txdbX,filenamesX,plotnameX,datadir1,datadir2){
  peakFiles<-lapply(filenamesX, function(filename) readPeakFile(paste0(datadir1,filename,'.bed')))
  names(peakFiles)<-filenamesX
  peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdbX,tssRegion=c(-3000, 3000), verbose=FALSE)
  tiff(paste0(datadir2,'peakAnno_',plotnameX,'.tiff'),height = 18, width = 15, units='cm', compression = "lzw", res = 300)
  print(plotAnnoBar(peakAnnoList))
  dev.off()

  tiff(paste0(datadir2,'distToTSS_',plotnameX,'.tiff'),height = 18, width = 15, units='cm', compression = "lzw", res = 300)
  print(plotDistToTSS(peakAnnoList))
  dev.off()
}

txdb1 <- TxDb.Hsapiens.UCSC.hg19.knownGene
filenames<-c("down_all3","down_RA_PBRA","down_PB_PBRA","down_PBRA","down_RA","down_PB","down_PB_RA","ns_all3","PBdown_RAup_PBRAdown","PBdown_RAup_PBRAns","PBdown_RAup_PBRAup","PBup_RAdown_PBRAdown","PBup_RAdown_PBRAns","PBup_RAdown_PBRAup","up_PB_RA","up_PB","up_RA","up_PBRA","up_PB_PBRA","up_RA_PBRA","up_all3")

plot_ChIPseeker_multiple(txdb1,filenames,'ChIPseeker',indir,outdir)


                    
