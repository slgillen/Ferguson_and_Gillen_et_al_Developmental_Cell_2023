#this script uses ChIPseeker to annotate the location of the H3K27ac broad peaks and to assign the closest gene to the peak for gene ontology analysis

# load libraries ----------------------------------------------------------
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)


# plot peak annotations with ChIPseeker -----------------------------------

plot_ChIPseeker_multiple<-function(txdbX,filenames,plotname,datadir){
  peakFiles<-lapply(filenames, function(filename) readPeakFile(paste0(datadir,filename,'.txt')))
  names(peakFiles)<-filenames
  peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb=txdbX,tssRegion=c(-3000, 3000), verbose=FALSE)
  tiff(paste0(datadir,'peakAnno_',plotname,'.tiff'),height = 5, width = 16.5, units='cm', compression = "lzw", res = 300)
  print(plotAnnoBar(peakAnnoList))
  dev.off()
  
  tiff(paste0(datadir,'distToTSS_',plotname,'.tiff'),height = 5, width = 16.5, units='cm', compression = "lzw", res = 300)
  print(plotDistToTSS(peakAnnoList))
  dev.off()
}

txdb1 <- TxDb.Hsapiens.UCSC.hg19.knownGene

BE2Cfiles<-c('BE2C_7dPBvControl_UP_FDR005Fold05','BE2C_7dPBvControl_ns_FDR05','BE2C_7dPBvControl_DOWN_FDR005Fold05')
IMR32files<-c('IMR32_5dPBvControl_UP_FDR005Fold05','IMR32_5dPBvControl_ns_FDR05','IMR32_5dPBvControl_DOWN_FDR005Fold05')
SY5Yfiles<-c('SY5Y_5dPBvControl_UP_FDR005Fold05','SY5Y_5dPBvControl_ns_FDR05','SY5Y_5dPBvControl_DOWN_FDR005Fold05')

filesets<-list(BE2Cfiles,IMR32files,SY5Yfiles)
plotnames<-c('BE2C_FC05','IMR32_FC05','SY5Y_FC05')

for(x in 1:length(filesets)){
  plot_ChIPseeker_multiple(txdb1,filesets[[x]],plotnames[x],'PB_H3K27ac_ChIPseq/data_output/')
}



# get peak annotation tables ----------------------------------------------

datanames<-c(BE2Cfiles,IMR32files,SY5Yfiles)

txdb1 <- TxDb.Hsapiens.UCSC.hg19.knownGene

datadir<-'PB_H3K27ac_ChIPseq/data_output/'

for(dataname in datanames){
  peak <- readPeakFile(paste0(datadir,filename,'.txt'))
  peakAnno.txdb <- annotatePeak(peak, tssRegion=c(-3000, 3000),TxDb=txdb1, annoDb="org.Hs.eg.db")
  peakAnno.txdb2<-data.frame(peakAnno.txdb,stringsAsFactors=FALSE)
  peakAnno.txdb2<-apply(peakAnno.txdb2,2,function(x) gsub('\t','',x))
  peakAnno.txdb2$annotation<-gsub('\t','',peakAnno.txdb2$annotation)
  write.table(peakAnno.txdb2,paste0(datadir,'peakAnno_',dataname,'.txt'),quote=FALSE,col.names=TRUE,row.names=FALSE,sep='\t')
  rm(peak,txdb1,peakAnno.txdb,peakAnno.txdb2)
}




# gene ontology analysis on proximal genes --------------------------------

# get proximal gene lists
# with filtering for a max distance between the H3K27ac broad peak and the gene TSS of 100kb
BE2Cup<-read.delim('peakAnno_BE2C_7dPBvControl_UP_FDR005Fold05.txt',stringsAsFactors = FALSE)
BE2Cup$distanceToTSS<-abs(BE2Cup$distanceToTSS)
BE2Cup_proxgenes<-unique(subset(BE2Cup,distanceToTSS<100000)$SYMBOL)

BE2Cdown<-read.delim('peakAnno_BE2C_7dPBvControl_DOWN_FDR005Fold05.txt',stringsAsFactors = FALSE)
BE2Cdown$distanceToTSS<-abs(BE2Cdown$distanceToTSS)
BE2Cdown_proxgenes<-unique(subset(BE2Cdown,distanceToTSS<100000)$SYMBOL)

IMR32up<-read.delim('peakAnno_IMR32_5dPBvControl_UP_FDR005Fold05.txt',stringsAsFactors = FALSE)
IMR32up$distanceToTSS<-abs(IMR32up$distanceToTSS)
IMR32up_proxgenes<-unique(subset(IMR32up,distanceToTSS<100000)$SYMBOL)

IMR32down<-read.delim('peakAnno_IMR32_5dPBvControl_DOWN_FDR005Fold05.txt',stringsAsFactors = FALSE)
IMR32down$distanceToTSS<-abs(IMR32down$distanceToTSS)
IMR32down_proxgenes<-unique(subset(IMR32down,distanceToTSS<100000)$SYMBOL)

SHSY5Yup<-read.delim('peakAnno_SHSY5Y_5dPBvControl_UP_FDR005Fold05.txt',stringsAsFactors = FALSE)
SHSY5Yup$distanceToTSS<-abs(SHSY5Yup$distanceToTSS)
SHSY5Yup_proxgenes<-unique(subset(SHSY5Yup,distanceToTSS<100000)$SYMBOL)

SHSY5Ydown<-read.delim('peakAnno_SHSY5Y_5dPBvControl_DOWN_FDR005Fold05.txt',stringsAsFactors = FALSE)
SHSY5Ydown$distanceToTSS<-abs(SHSY5Ydown$distanceToTSS)
SHSY5Ydown_proxgenes<-unique(subset(SHSY5Ydown,distanceToTSS<100000)$SYMBOL)


#get set of genes that are proximal to increased H3K27ac deposition across cell lines
overlap_up<-Reduce(intersect, list(BE2Cup_proxgenes,IMR32up_proxgenes,SY5Yup_proxgenes))

#get set of genes that are proximal to decreased H3K27ac deposition across cell lines
overlap_down<-Reduce(intersect, list(BE2Cdown_proxgenes,IMR32down_proxgenes,SY5Ydown_proxgenes))

# run gene ontology analysis
genelists<-list(overlap_up,overlap_down)
setnames<-list('overlap_up','overlap_down')

GOtypes<-c('BP')
for(GOtype in GOtypes){
  for(i in 1:length(genelists)){
    geneList<-genelists[[i]]
    ego <- enrichGO(geneList, OrgDb = "org.Hs.eg.db", ont=GOtype, pAdjustMethod='fdr', pvalueCutoff=0.05,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
    ego_df<-data.frame(ego,stringsAsFactors = FALSE)
    ego_df<-ego_df[order(ego_df$p.adjust),]
    write.csv(ego_df[,c(2,3,4,6,8)],file = paste0('PB_H3K27ac_ChIPseq/output_data/',setnames[i],'_',GOtype,'.csv'),quote=FALSE,row.names=FALSE)
    ego_df$neglog10padj<--log(ego_df$p.adjust,10)
    ego_df<-ego_df[order(ego_df$p.adjust,decreasing=TRUE),] #order table by p-value
    
    ego_df$Description<-factor(ego_df$Description,level=ego_df$Description)
    ego_df$gene_ratio<-sapply(ego_df$GeneRatio, function(x) as.numeric(strsplit(x,split='/')[[1]][1])/as.numeric(strsplit(x,split='/')[[1]][2]))
    
    #if more than 10 significant terms, then plot the 10 most significant
    if(nrow(ego_df)>=10){
      png(paste0('PB_H3K27ac_ChIPseq/output_data/',setnames[i],'_',GOtype,'.png'),width=15,height=5,units = 'cm',res=200)
      print(ggplot(ego_df[((nrow(ego_df)-9):nrow(ego_df)),],aes(x=neglog10padj,y=Description,fill=gene_ratio))+geom_bar(stat='identity')+theme_classic()+
              scale_x_continuous(expand=c(0,0))+labs(fill = 'Gene Ratio')+xlab('-log10 adjusted p-val')+ylab('')+theme(axis.text=element_text(size=12),axis.title=element_text(size=12)))
      dev.off()
    }else{
      png(paste0('PB_H3K27ac_ChIPseq/output_data/',setnames[i],'_',GOtype,'.png'),width=30,height=10,units = 'cm',res=200)
      print(ggplot(ego_df,aes(x=neglog10padj,y=Description,fill=gene_ratio))+geom_bar(stat='identity')+theme_classic()+
              scale_x_continuous(expand=c(0,0))+labs(fill = 'Gene Ratio')+xlab('-log10 adjusted p-val')+ylab('')+theme(axis.text=element_text(size=12),axis.title=element_text(size=12)))
      dev.off()
    }
    rm(ego,ego_df,geneList)
  }
}

