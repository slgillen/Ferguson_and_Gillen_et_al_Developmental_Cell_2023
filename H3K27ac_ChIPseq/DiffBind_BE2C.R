#this scipt conducts DiffBind analysis to identify regions with significantly altered H3K27ac deposition between PB-treated and control samples in the SK-N-BE(2)C data

# load libraries ----------------------------------------------------------
library(DiffBind)
library(BiocParallel)


# read in data ------------------------------------------------------------
#broad peaks were called with MACS2 and input used as a reference
samples <- read.csv("PB_H3K27ac_ChIPseq/data/BE2C_DiffBind_SampleSheet_vInput.csv",stringsAsFactors=FALSE)
print(names(samples))

BE2C_data <- dba(sampleSheet=samples,minOverlap=1)


# black list removal ------------------------------------------------------
BE2C_data<-dba.blacklist(BE2C_data,blacklist=DBA_BLACKLIST_HG19,cores=10)


# check sample correlations -----------------------------------------------
pdf(paste0(datadir,'PB_H3K27ac_ChIPseq/output_figures/BE2C_correlations_ConsensusNoSummits_vInput_broad.pdf'),width=500,height=400)
plot(BE2C_data)
dev.off()


# get consensus peak set across samples and conditions --------------------
info<-dba.show(BE2C_data) 
print(info)

# check peak overlaps across replicates
olap.rate.control <- dba.overlap(BE2C_data,BE2C_data$masks$control,mode=DBA_OLAP_RATE) 
print(olap.rate.control)
olap.rate.PB <- dba.overlap(BE2C_data,BE2C_data$masks$PB5d,mode=DBA_OLAP_RATE) 
print(olap.rate.PB)

# check peak overlaps across all samples
olap.rate.all <- dba.overlap(BE2C_data,mode=DBA_OLAP_RATE) 
print(olap.rate.all) 

# get consensus peak sets for each condition
#keep if in at least 3 of the 5 control replicates
control_consensus_1<-dba.peakset(BE2C_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_data$masks$control, minOverlap=3) 
control_consensus_1x<-control_consensus_1[,1:3]
nrow(control_consensus_1x) 
write.table(control_consensus_1x,paste0(datadir,'PB_H3K27ac_ChIPseq/output_data/BE2C_control_ConsensusNoSummits_vInput_broad_in3of5.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PB_consensus_1<-dba.peakset(BE2C_data, bRetrieve=T, DataType=DBA_DATA_FRAME,BE2C_data$masks$PB, minOverlap=3) 
PB_consensus_1x<-PB_consensus_1[,1:3]
nrow(PB_consensus_1x) 
write.table(PB_consensus_1x,paste0(datadir,'PB_H3K27ac_ChIPseq/output_data/BE2C_PB_ConsensusNoSummits_vInput_broad_in3of5.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


# get full consensus peakset across conditions
BE2C_consensus <- dba.peakset(BE2C_data, consensus = DBA_CONDITION, minOverlap=3)
ExptConsensus <- dba(BE2C_consensus, mask=BE2C_consensus$masks$Consensus,minOverlap=1)
ExptConsensus
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE,DataType=DBA_DATA_FRAME,minOverlap=1)
ConsensusPeaksx<-ConsensusPeaks[,1:3]
nrow(ConsensusPeaksx) 
write.table(ConsensusPeaksx,paste0(datadir,'PB_H3K27ac_ChIPseq/output_data/BE2C_ConsensusNoSummits_vInput_broad_bothin3of5.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


# get peakset counts ------------------------------------------------------
register(MulticoreParam(4))
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE)
BE2C_data_count <- dba.count(BE2C_data,peaks=ConsensusPeaks,bParallel=TRUE,summits=FALSE) 


# plot PCA ----------------------------------------------------------------
tiff(file='PB_H3K27ac_ChIPseq/output_figures/BE2C_PCA_normalised_ConsensusNoSummits_vInput_broad_nolabel.tiff',width=1400,height=1300,res=300)
dba.plotPCA(BE2C_data_count_v2,DBA_CONDITION,vColors=c('grey48','red3'),score=DBA_SCORE_NORMALIZED)
dev.off()


# get FRiP details; FRiP = fraction of reads in the peaks -----------------
info <- dba.show(BE2C_data_count_v2)
write.table(data.frame(info,stringsAsFactors=FALSE),'PB_H3K27ac_ChIPseq/output_figures/info_BE2C_vInput_broad.txt',col.names=TRUE,row.names=FALSE,quote=FALSE)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
print(libsizes)


# normalisation -----------------------------------------------------------
BE2C_data_normalize<-dba.normalize(BE2C_data_count_v2,method=DBA_DESEQ2) 


# conduct differential analysis -------------------------------------------

# sort contrast
BE2C_data_contrast <- dba.contrast(BE2C_data_normalize,design="~Replicate+Condition") 

# run analysis
BE2C_data_analyze <- dba.analyze(BE2C_data_contrast,method = DBA_DESEQ2,bGreylist=FALSE) 
dba.show(BE2C_data_analyze, bContrasts=TRUE)

# check differential sites
BE2C_data_analyze_DiffSites <- dba.report(BE2C_data_analyze,method = DBA_DESEQ2) 
print(sum(BE2C_data_analyze_DiffSites$Fold>0 & BE2C_data_analyze_DiffSites$FDR<0.05))
print(sum(BE2C_data_analyze_DiffSites$Fold<0 & BE2C_data_analyze_DiffSites$FDR<0.05))



# write outputs -----------------------------------------------------------

# DiffBind output table
report <- dba.report(BE2C_data_analyze, contrast = 1, th = 1, bFlip = F)
write.table(report, paste0("PB_H3K27ac_ChIPseq/output_data/BE2C_7dPBvControl_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

# plot data as volcano plot
png('PB_H3K27ac_ChIPseq/output_figures/BE2C_volcano_afteranalysis_ConsensusNoSummits_vInput_broad.png',width=500,height=400)
dba.plotVolcano(BE2C_data_analyze) 
dev.off()

report_df<-data.frame(report)
nrow(report_df)
nrow(subset(report_df,FDR<0.05 & Fold>0.5))
nrow(subset(report_df,FDR<0.05 & Fold<(-0.5)))
nrow(subset(report_df,FDR>0.5))



# get sub-group files for downstream plotting ---------------

# for input into ChIPseeker
write.table(subset(report_df,FDR<0.05 & Fold>0.5), paste0("BE2C_5dPBvControl_UP_FDR005Fold05.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(subset(report_df,FDR<0.05 & Fold<(-0.5)), paste0("BE2C_5dPBvControl_DOWN_FDR005Fold05.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(subset(report_df,FDR>0.5), paste0("BE2C_5dPBvControl_ns_FDR05.txt"), sep = "\t", col.names = T, row.names = F, quote = F)


# for input into plotHeatmap
BE2C_broad<-report_df
BE2C_broad$name<-apply(BE2C_broad[,1:3],1,function(x) paste(x[1:3],collapse=':'))
BE2C_broad$name<-gsub(' ','',BE2C_broad$name)

#order outputs by FDR
BE2C_broad<-BE2C_broad[order(BE2C_broad$FDR,decreasing=FALSE),]
BE2C_broad_sigup<-subset(BE2C_broad,Fold>0.5 & FDR<0.05)
nrow(BE2C_broad_sigup)
BE2C_broad_sigdown<-subset(BE2C_broad,Fold<(-0.5) & FDR<0.05)
nrow(BE2C_broad_sigdown)
BE2C_broad_ns<-subset(BE2C_broad,FDR>0.5)
nrow(BE2C_broad_ns)

write.table(BE2C_broad_sigup[,c(1,2,3,4,5,11)],'BE2C_broadpeak_up05.bed',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(BE2C_broad_sigdown[,c(1,2,3,4,5,11)],'BE2C_broadpeak_down05.bed',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(BE2C_broad_ns[,c(1,2,3,4,5,11)],'BE2C_broadpeak_ns05.bed',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)


