#this scipt conducts DiffBind analysis to identify regions with significantly altered H3K27ac deposition between DMSO/RA/PB/PB+RA treated samples in the SK-N-BE(2)C data

#using R version: 4.1.2

#load libraries
library(DiffBind) #version 3.4
library(BiocParallel)


# read in data ------------------------------------------------------------
samples <- read.csv(paste0(datadir,"PBandRA_H3K27ac_ChIPseq/data/PBandRA_DiffBind_SampleSheet_vInput.csv"),stringsAsFactors=FALSE) #narrow just means MACS2
print(names(samples))

BE2C_PBandRA_data <- dba(sampleSheet=samples,minOverlap=1)

# blacklist filtering ------------------------------------------------------------
BE2C_PBandRA_data<-dba.blacklist(BE2C_PBandRA_data,blacklist=DBA_BLACKLIST_HG19,cores=12)


# check sample correlations ------------------------------------------------------------
tiff(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_correlations.tiff'),width=1400,height=1400,res=300)
plot(BE2C_PBandRA_data)
dev.off() 


# get peak overlap rates across conditions --------------------
info<-dba.show(BE2C_PBandRA_data) 
print(info)

olap.rate.DMSO <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$DMSO,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap DMSO',file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
write(olap.rate.DMSO,file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
print(olap.rate.DMSO)

olap.rate.RA <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$RA,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap RA',file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
write(olap.rate.RA,file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
print(olap.rate.RA)

olap.rate.PB <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$PB,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap PB',file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
write(olap.rate.PB,file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
print(olap.rate.PB)

olap.rate.PBRA <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$PBRA,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap PBRA',file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
write(olap.rate.PBRA,file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
print(olap.rate.PBRA)

olap.rate.all <- dba.overlap(BE2C_PBandRA_data,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap all',file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
write(olap.rate.all,file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/overlap_rates.txt'),append=TRUE)
print(olap.rate.all)


# get consensus peak set across samples and conditions --------------------

# keep if in at least 2 of the 4 replicates in a condition
DMSO_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$DMSO, minOverlap=2) 
DMSO_consensus_1x<-DMSO_consensus_1[,1:3]
nrow(DMSO_consensus_1x) 
write.table(DMSO_consensus_1x,paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/DMSO_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

RA_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$RA, minOverlap=2) 
RA_consensus_1x<-RA_consensus_1[,1:3]
nrow(RA_consensus_1x) 
write.table(RA_consensus_1x,paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/RA_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PB_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$PB, minOverlap=2) 
PB_consensus_1x<-PB_consensus_1[,1:3]
nrow(PB_consensus_1x) 
write.table(PB_consensus_1x,paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PB_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PBRA_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$PBRA, minOverlap=2) 
PBRA_consensus_1x<-PBRA_consensus_1[,1:3]
nrow(PBRA_consensus_1x) 
write.table(PBRA_consensus_1x,paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRA_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


# get full consensus peakset across conditions
BE2C_PBandRA_consensus <-  dba.peakset(BE2C_PBandRA_data, consensus = DBA_CONDITION, minOverlap=2)
ExptConsensus <-  dba(BE2C_PBandRA_consensus, mask=BE2C_PBandRA_consensus$masks$Consensus,minOverlap=1) 
ExptConsensus
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE,DataType=DBA_DATA_FRAME,minOverlap=1)
ConsensusPeaksx<-ConsensusPeaks[,1:3]
nrow(ConsensusPeaksx) 
write.table(ConsensusPeaksx,paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_ConsensusPeaks_allin2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


# plots with consensus data --------------------
tiff(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_overlapVenn_ConsensusPeaks.tiff'),width=1400,height=1400,res=300)
dba.plotVenn(ExptConsensus,ExptConsensus$masks$Consensus)
dev.off()

tiff(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_ConsensusPeaks.tiff'),width=1400,height=1400,res=300)
dba.plotPCA(BE2C_PBandRA_consensus,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'))
dev.off() 


# get peakset counts for consensus regions defined --------------------

#sort config
BE2C_PBandRA_data$config$singleEnd <- FALSE #using paired-end data
BE2C_PBandRA_data$config$cores <- 6
BE2C_PBandRA_data$config

#get the counts
ConsensusPeaks_v2 <- dba.peakset(ExptConsensus, bRetrieve=TRUE)
BE2C_PBandRA_data_count <- dba.count(BE2C_PBandRA_data,peaks=ConsensusPeaks_v2,bParallel=TRUE,summits=FALSE,filter=1,minCount=0,bUseSummarizeOverlaps=TRUE) 
BE2C_PBandRA_data_count

reads <- dba.peakset(BE2C_PBandRA_data_count, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
apply(reads[,4:19],2,function(x) sum(x))

#write output
write.csv(reads, file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_inconsensuspeaks_reads.csv'))


# PCA plots with consensus data ----------------------------------------------

#FRiP = fraction of reads in the peaks
info <- dba.show(BE2C_PBandRA_data_count)
write.table(data.frame(info,stringsAsFactors=FALSE),paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/info_BE2C_PBandRA.txt'),col.names=TRUE,row.names=FALSE,quote=FALSE)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
print(libsizes)

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_normalised_Consensus.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_normalised_Consensus.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()


#check PB v DMSO
dba_ss <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$DMSO)
dba_ss
dba_ss2 <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$PB)
dba_ss2
myDBA <- dba.peakset(dba_ss,peaks=dba_ss2)
myDBA

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_normalised_PBvDMSO.png'),width=375,height=320)
dba.plotPCA(myDBA,components=1:2,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','red3'),score=DBA_SCORE_NORMALIZED)
dev.off()

#check PB+RA v RA
dba_ss <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$RA)
dba_ss
dba_ss2 <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$PBRA)
dba_ss2
myDBA <- dba.peakset(dba_ss,peaks=dba_ss2)
myDBA

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_normalised_PBRAvRA.png'),width=375,height=320)
dba.plotPCA(myDBA,components=1:2,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('royalblue3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_normalised_Consensus_PC34.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,components=3:4,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/BE2C_PBandRA_PCA_normalised_Consensus_v2.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,DBA_CONDITION,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

      
# data normalisation ----------------------------------------------

BE2C_PBandRA_normalize<-dba.normalize(BE2C_PBandRA_data_count,method=DBA_DESEQ2)
print(BE2C_PBandRA_normalize$norm)
normlibs <- cbind(FullLibSize=BE2C_PBandRA_normalize$norm$DESeq2$lib.sizes, NormFacs=BE2C_PBandRA_normalize$norm$DESeq2$norm.facs,NormLibSize=round(BE2C_PBandRA_normalize$norm$DESeq2$lib.sizes/BE2C_PBandRA_normalize$norm$DESeq2$norm.facs))
rownames(normlibs) <- info$ID
print(normlibs)


# conduct differential analysis ----------------------------------------------
BE2C_PBandRA_contrast <- dba.contrast(BE2C_PBandRA_normalize,design="~Replicate+Condition")

BE2C_PBandRA_analyze <- dba.analyze(BE2C_PBandRA_contrast,method = DBA_DESEQ2,bGreylist=FALSE,bBlacklist=FALSE) 
dba.show(BE2C_PBandRA_analyze, bContrasts=TRUE)
write.table(dba.show(BE2C_PBandRA_analyze, bContrasts=TRUE),file=paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBandRA_analyze_info.txt'),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)




# write outputs for the different comparisons ----------------------------------------------
#bFLIP=T or bFLIP=F to sort comparisons direction

report <- dba.report(BE2C_PBandRA_analyze, contrast = 7, th = 1, bFlip = T)
write.table(report, paste0("PBandRA_H3K27ac_ChIPseq/DiffBind/RAvDMSO_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 8, th = 1, bFlip = T)
write.table(report, paste0("PBandRA_H3K27ac_ChIPseq/DiffBind/PBvDMSO_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 9, th = 1, bFlip = T)
write.table(report, paste0("PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvDMSO_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 10, th = 1, bFlip = T)
write.table(report, paste0("PBandRA_H3K27ac_ChIPseq/DiffBind/PBvRA_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 11, th = 1, bFlip = T)
write.table(report, paste0("PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvRA_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 12, th = 1, bFlip = F) 
write.table(report, paste0("PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvPB_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

      
# visualisations of the differential ----------------------------------------------
#should make own volcano plots

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/RAvDMSO_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, contrast=7, label=DBA_REPLICATE,vColors=c('grey48','royalblue3'))
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/RAvDMSO_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=7)
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/RAvDMSO_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, bFlip = T, contrast=7)
dev.off()

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBvDMSO_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=8, label=DBA_REPLICATE,vColors=c('grey48','red3'))
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBvDMSO_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=8)
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBvDMSO_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=8)
dev.off()

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvDMSO_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=9, label=DBA_REPLICATE,vColors=c('grey48','purple'))
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvDMSO_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=9)
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvDMSO_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=9)
dev.off()


png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBvRA_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=10, label=DBA_REPLICATE,vColors=c('royalblue3','red3'))
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBvRA_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=10)
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBvRA_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=10)
dev.off()


png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvRA_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=11, label=DBA_REPLICATE,vColors=c('royalblue3','purple'))
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvRA_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=11)
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvRA_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=11) 
dev.off()

png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvPB_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, contrast=12, label=DBA_REPLICATE,vColors=c('purple','red3'))
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvPB_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, contrast=12)
dev.off()
png(paste0('PBandRA_H3K27ac_ChIPseq/DiffBind/PBRAvPB_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=12) 
dev.off()






##################### other data normalisation PCA #####################

print(BE2C_PBandRA_data_count$SN) #can update this using more accurate rounding from the raw reads data above if required - only necessary with TF data that has much lower FRiP scores?


#normalise the data and replot PCA - with DESeq2 norm #i.e. RLE norm?
BE2C_PBandRA_data_DESEQnorm<-dba.normalize(BE2C_PBandRA_data_count,method=DBA_DESEQ2,minCount=1,peaks=ConsensusPeaks_v2) 

tiff(paste0(datadir,'PBandRA_data_DESEQnorm_PCA.tiff'),width=1400,height=1400,res=300)
dba.plotPCA(BE2C_PBandRA_data_DESEQnorm,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'))
dev.off()

BE2C_PBandRA_data_FRIPnorm<-dba.normalize(BE2C_PBandRA_data_count,normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_PEAKREADS) 

tiff(paste0(datadir,'PBandRA_FRiPnorm_PCA.tiff'),width=1400,height=1400,res=300)
dba.plotPCA(BE2C_PBandRA_data_FRIPnorm,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'))
dev.off()


