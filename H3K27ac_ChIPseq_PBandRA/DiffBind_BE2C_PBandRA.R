#using R version: 4.1.2

#load libraries
library(DiffBind) #version 3.4
library(BiocParallel)


# read in data ------------------------------------------------------------
samples <- read.csv(paste0(datadir,"PB_H3K27ac_ChIPseq/data/PBandRA_DiffBind_SampleSheet_vInput.csv"),stringsAsFactors=FALSE) #narrow just means MACS2
print(names(samples))

BE2C_PBandRA_data <- dba(sampleSheet=samples,minOverlap=1)

# blacklist filtering ------------------------------------------------------------
BE2C_PBandRA_data<-dba.blacklist(BE2C_PBandRA_data,blacklist=DBA_BLACKLIST_HG19,cores=12)


##################### plot correlations <- should try change settings on this for value range at least #####################
tiff(paste0(datadir,'BE2C_PBandRA_correlations.tiff'),width=1400,height=1400,res=300)
plot(BE2C_PBandRA_data)
dev.off() #should retry correlations at later stage on normalised within peak count data - this could just show calling diffs

info<-dba.show(BE2C_PBandRA_data) 
print(info)

tiff(paste0(datadir,'BE2C_PBandRA_PCA.tiff'),width=1400,height=1400,res=300)
dba.plotPCA(BE2C_PBandRA_data,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'))
dev.off() #should retry PCA at later stage on normalised within peak count data - this could just show calling diffs


##################### get overlap rates #####################
olap.rate.DMSO <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$DMSO,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap DMSO',file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
write(olap.rate.DMSO,file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
print(olap.rate.DMSO)

olap.rate.RA <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$RA,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap RA',file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
write(olap.rate.RA,file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
print(olap.rate.RA)

olap.rate.PB <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$PB,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap PB',file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
write(olap.rate.PB,file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
print(olap.rate.PB)

olap.rate.PBRA <- dba.overlap(BE2C_PBandRA_data,BE2C_PBandRA_data$masks$PBRA,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap PBRA',file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
write(olap.rate.PBRA,file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
print(olap.rate.PBRA)

olap.rate.all <- dba.overlap(BE2C_PBandRA_data,mode=DBA_OLAP_RATE) #overlap rate of peaks - number that occur in 1,2.. etc samples
write('overlap all',file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
write(olap.rate.all,file=paste0(datadir,'overlap_rates.txt'),append=TRUE)
print(olap.rate.all) #diffbind is retaining for analysis that is in 2 replicates of everything...
#actually want to keep everything within control or PB consensus - maybe limit these to 2 instead?





##################### consensus peak definitions #####################
#default is in at least 2 of all the samples

#keep if in at least 3 of the 4 DMSO replicates
DMSO_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$DMSO, minOverlap=3) 
DMSO_consensus_1x<-DMSO_consensus_1[,1:3]
nrow(DMSO_consensus_1x) 
write.table(DMSO_consensus_1x,paste0(datadir,'DMSO_peaks_in3of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

RA_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$RA, minOverlap=3) 
RA_consensus_1x<-RA_consensus_1[,1:3]
nrow(RA_consensus_1x) 
write.table(RA_consensus_1x,paste0(datadir,'RA_peaks_in3of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PB_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$PB, minOverlap=3) 
PB_consensus_1x<-PB_consensus_1[,1:3]
nrow(PB_consensus_1x) 
write.table(PB_consensus_1x,paste0(datadir,'PB_peaks_in3of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PBRA_consensus_1<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$PBRA, minOverlap=3) 
PBRA_consensus_1x<-PBRA_consensus_1[,1:3]
nrow(PBRA_consensus_1x) 
write.table(PBRA_consensus_1x,paste0(datadir,'PBRA_peaks_in3of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")



#also get less stringent keep in at least 2 of the 4 control replicates in case needed
DMSO_consensus_2<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$DMSO, minOverlap=2) 
DMSO_consensus_2x<-DMSO_consensus_2[,1:3]
nrow(DMSO_consensus_2x) 
write.table(DMSO_consensus_2x,paste0(datadir,'DMSO_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

RA_consensus_2<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$RA, minOverlap=2) 
RA_consensus_2x<-RA_consensus_2[,1:3]
nrow(RA_consensus_2x) 
write.table(RA_consensus_2x,paste0(datadir,'RA_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PB_consensus_2<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$PB, minOverlap=2) 
PB_consensus_2x<-PB_consensus_2[,1:3]
nrow(PB_consensus_2x) 
write.table(PB_consensus_2x,paste0(datadir,'PB_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PBRA_consensus_2<-dba.peakset(BE2C_PBandRA_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_PBandRA_data$masks$PBRA, minOverlap=2) 
PBRA_consensus_2x<-PBRA_consensus_2[,1:3]
nrow(PBRA_consensus_2x) 
write.table(PBRA_consensus_2x,paste0(datadir,'PBRA_peaks_in2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")





##################### get full consensus peak set - have as called in at least 2 of 4 replicate of any given condition #####################

###### ************* NOTE: check this stage now makes sense when using the multiple conditions ******** #######

BE2C_PBandRA_consensus <-  dba.peakset(BE2C_PBandRA_data, consensus = DBA_CONDITION, minOverlap=2)
ExptConsensus <-  dba(BE2C_PBandRA_consensus, mask=BE2C_PBandRA_consensus$masks$Consensus,minOverlap=1) 
ExptConsensus
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE,DataType=DBA_DATA_FRAME,minOverlap=1)
ConsensusPeaksx<-ConsensusPeaks[,1:3]
nrow(ConsensusPeaksx) 
write.table(ConsensusPeaksx,paste0(datadir,'BE2C_PBandRA_ConsensusPeaks_bothin2of4.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


##################### consensus plotting #####################
tiff(paste0(datadir,'BE2C_PBandRA_overlapVenn_ConsensusPeaks.tiff'),width=1400,height=1400,res=300)
dba.plotVenn(ExptConsensus,ExptConsensus$masks$Consensus)
dev.off()

tiff(paste0(datadir,'BE2C_PBandRA_PCA_ConsensusPeaks.tiff'),width=1400,height=1400,res=300)
dba.plotPCA(BE2C_PBandRA_consensus,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'))
dev.off() 


##################### get read counts (compare PE and SE settings? - does SE lead to better stats due to double count?) #####################

ConsensusPeaks_v2 <- dba.peakset(ExptConsensus, bRetrieve=TRUE)
BE2C_PBandRA_data$config$singleEnd <- FALSE #would be good to compare FRiP scores with and without this
BE2C_PBandRA_data$config$cores <- 6
BE2C_PBandRA_data$config
#don't set summits because small range excludes key aspects?
#now states reads will be counted as paired-end
#filter and mincount values changed
#NOTE: default of bRemoveDuplicates=FALSE
BE2C_PBandRA_data_count <- dba.count(BE2C_PBandRA_data,peaks=ConsensusPeaks_v2,bParallel=TRUE,summits=FALSE,filter=1,minCount=0,bUseSummarizeOverlaps=TRUE) #only works with summits=TRUE for some reason, not when setting specific size range
BE2C_PBandRA_data_count #copied this info to text file #FRiPs are more variable - lowest is PB R4

reads <- dba.peakset(BE2C_PBandRA_data_count, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
apply(reads[,4:19],2,function(x) sum(x))
write.csv(reads, file=paste0(datadir,"BE2C_inconsensuspeaks_reads.csv"))


# #testing comparison to SE setting
# BE2C_PBandRA_data_SE<-BE2C_PBandRA_data
# BE2C_PBandRA_data_SE$config$singleEnd <- TRUE #would be good to compare FRiP scores with and without this
# #don't set summits because small range excludes key aspects?
# #now states reads will be counted as paired-end
# BE2C_PBandRA_data_count_SE <- dba.count(BE2C_PBandRA_data_SE,peaks=ConsensusPeaks_v2,bParallel=TRUE,summits=FALSE,filter=0,minCount=1,bUseSummarizeOverlaps=TRUE) #only works with summits=TRUE for some reason, not when setting specific size range
# BE2C_PBandRA_data_count_SE #copied this info to text file

# reads_SE <- dba.peakset(BE2C_PBandRA_data_SE, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
# apply(reads_SE[,4:19],2,function(x) sum(x))
# write.csv(reads_SE, file=paste0(datadir,"BE2C_inconsensuspeaks_reads_SE.csv"))



###################### save the dba object ########################
savefile <- dba.save(BE2C_PBandRA_data_count,dir=datadir,'BE2C_PBandRA_data_count')
savefile
BE2C_PBandRA_data_countX <- dba.load(dir=datadir,'BE2C_PBandRA_data_count')
unlink(savefile)


##################### data normalisation and PCA plot #####################

#FRiP = fraction of reads in the peaks
info <- dba.show(BE2C_PBandRA_data_count)
write.table(data.frame(info,stringsAsFactors=FALSE),paste0(datadir,'info_BE2C_PBandRA_vInput_broad.txt'),col.names=TRUE,row.names=FALSE,quote=FALSE)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
print(libsizes)

png(paste0(datadir,'BE2C_PBandRA_PCA_normalised_ConsensusNoSummits_vInput_broad.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

png(paste0(datadir,'BE2C_PBandRA_PCA_normalised_ConsensusNoSummits_vInput_broad.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

dba_ss <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$DMSO)
dba_ss
dba_ss2 <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$PB)
dba_ss2
myDBA <- dba.peakset(dba_ss,peaks=dba_ss2)
myDBA

png(paste0(datadir,'BE2C_PBandRA_PCA_normalised_PBvDMSOcheck.png'),width=375,height=320)
dba.plotPCA(myDBA,components=1:2,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','red3'),score=DBA_SCORE_NORMALIZED)
dev.off()


dba_ss <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$RA)
dba_ss
dba_ss2 <- dba(BE2C_PBandRA_data_count, mask=BE2C_PBandRA_data_count$masks$PBRA)
dba_ss2
myDBA <- dba.peakset(dba_ss,peaks=dba_ss2)
myDBA

png(paste0(datadir,'BE2C_PBandRA_PCA_normalised_PBRAvRAcheck.png'),width=375,height=320)
dba.plotPCA(myDBA,components=1:2,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('royalblue3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

png(paste0(datadir,'BE2C_PBandRA_PCA_normalised_ConsensusNoSummits_vInput_broad_PC34.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,components=3:4,DBA_CONDITION,label=DBA_REPLICATE,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()

png(paste0(datadir,'BE2C_PBandRA_PCA_normalised_ConsensusNoSummits_vInput_broad_nolabel.png'),width=375,height=320)
dba.plotPCA(BE2C_PBandRA_data_count,DBA_CONDITION,vColors=c('grey48','royalblue3','red3','purple'),score=DBA_SCORE_NORMALIZED)
dev.off()



#default normalisation is sequencing depth <- can use spike-ins at this point
#could normalise at this point based on FRiP instead??? - have done for TFs
BE2C_PBandRA_normalize<-dba.normalize(BE2C_PBandRA_data_count,method=DBA_DESEQ2) #setting bRetrieve=TRUE leads to having no output
#to see DESeq2 norm factors
print(BE2C_PBandRA_normalize$norm)
#could use DESeq2 generated normalisation factors
normlibs <- cbind(FullLibSize=BE2C_PBandRA_normalize$norm$DESeq2$lib.sizes, NormFacs=BE2C_PBandRA_normalize$norm$DESeq2$norm.facs,NormLibSize=round(BE2C_PBandRA_normalize$norm$DESeq2$lib.sizes/BE2C_PBandRA_normalize$norm$DESeq2$norm.facs))
rownames(normlibs) <- info$ID
print(normlibs)


##################### analyse - consider all different condition comparisons #####################

#sort contrast / design
#how best to do for comparisons across conditions? - filter down to just the comparisons each time? What are the output formats via DiffBind?
BE2C_PBandRA_contrast <- dba.contrast(BE2C_PBandRA_normalize,design="~Replicate+Condition") #Replicate as blocking factor


#check if differentially bound sites reflect same sites as narrowPeak files
BE2C_PBandRA_analyze <- dba.analyze(BE2C_PBandRA_contrast,method = DBA_DESEQ2,bGreylist=FALSE,bBlacklist=FALSE) #already done at earlier stage... why is this still applying black/grey lists? 
dba.show(BE2C_PBandRA_analyze, bContrasts=TRUE)
write.table(dba.show(BE2C_PBandRA_analyze, bContrasts=TRUE),file=paste0(datadir,'PBandRA_analyze_info.txt'),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
#NOTE: contrast is Group v Group2



##################### get differential #####################
#bFLIP=T or bFLIP=F to sort comparisons

report <- dba.report(BE2C_PBandRA_analyze, contrast = 7, th = 1, bFlip = T)
write.table(report, paste0("RAvDMSO_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 8, th = 1, bFlip = T)
write.table(report, paste0("PBvDMSO_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 9, th = 1, bFlip = T)
write.table(report, paste0("PBRAvDMSO_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 10, th = 1, bFlip = T)
write.table(report, paste0("PBvRA_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 11, th = 1, bFlip = T)
write.table(report, paste0("PBRAvRA_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

report <- dba.report(BE2C_PBandRA_analyze, contrast = 12, th = 1, bFlip = F) 
write.table(report, paste0("PBRAvPB_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
rm(report)

##################### visualise differential #####################
#should make own volcano plots
#how to sort the flip for these?? use bFLIP=TRUE/FALSE here also

png(paste0(datadir,'RAvDMSO_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, contrast=7, label=DBA_REPLICATE,vColors=c('grey48','royalblue3'))
dev.off()
png(paste0(datadir,'RAvDMSO_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=7)
dev.off()
png(paste0(datadir,'RAvDMSO_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, bFlip = T, contrast=7)
dev.off()

png(paste0(datadir,'PBvDMSO_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=8, label=DBA_REPLICATE,vColors=c('grey48','red3'))
dev.off()
png(paste0(datadir,'PBvDMSO_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=8)
dev.off()
png(paste0(datadir,'PBvDMSO_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=8)
dev.off()

png(paste0(datadir,'PBRAvDMSO_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=9, label=DBA_REPLICATE,vColors=c('grey48','purple'))
dev.off()
png(paste0(datadir,'PBRAvDMSO_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=9)
dev.off()
png(paste0(datadir,'PBRAvDMSO_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=9)
dev.off()


png(paste0(datadir,'PBvRA_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=10, label=DBA_REPLICATE,vColors=c('royalblue3','red3'))
dev.off()
png(paste0(datadir,'PBvRA_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=10)
dev.off()
png(paste0(datadir,'PBvRA_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=10)
dev.off()


png(paste0(datadir,'PBRAvRA_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, bFlip = T, contrast=11, label=DBA_REPLICATE,vColors=c('royalblue3','purple'))
dev.off()
png(paste0(datadir,'PBRAvRA_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, bFlip = T, contrast=11)
dev.off()
png(paste0(datadir,'PBRAvRA_volcano_afteranalysis.png'),width=500,height=400)
dba.plotVolcano(BE2C_PBandRA_analyze, contrast=11) 
dev.off()

png(paste0(datadir,'PBRAvPB_PCA_afteranalysis.png'),width=500,height=400)
dba.plotPCA(BE2C_PBandRA_analyze, contrast=12, label=DBA_REPLICATE,vColors=c('purple','red3'))
dev.off()
png(paste0(datadir,'PBRAvPB_MAplot_afteranalysis.png'),width=500,height=400)
dba.plotMA(BE2C_PBandRA_analyze, contrast=12)
dev.off()
png(paste0(datadir,'PBRAvPB_volcano_afteranalysis.png'),width=500,height=400)
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


