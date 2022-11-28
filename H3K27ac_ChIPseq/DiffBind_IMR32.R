#this scipt conducts DiffBind analysis to identify regions with significantly altered H3K27ac deposition between PB-treated and control samples in the IMR-32 data

# load libraries ----------------------------------------------------------
library(DiffBind)
library(BiocParallel)


# read in data ------------------------------------------------------------
#broad peaks were called with MACS2 and input used as a reference
samples <- read.csv("PB_H3K27ac_ChIPseq/data/IMR32_DiffBind_SampleSheet_vInput.csv",stringsAsFactors=FALSE)
print(names(samples))

IMR32_data <- dba(sampleSheet=samples,minOverlap=1)


# black list removal ------------------------------------------------------
IMR32_data<-dba.blacklist(IMR32_data,blacklist=DBA_BLACKLIST_HG19,cores=10)


# check sample correlations -----------------------------------------------
pdf(paste0(datadir,'PB_H3K27ac_ChIPseq/output_figures/IMR32_correlations_ConsensusNoSummits_vInput_broad.pdf'),width=500,height=400)
plot(IMR32_data)
dev.off()


# get consensus peak set across samples and conditions --------------------
info<-dba.show(IMR32_data) 
print(info)

# check peak overlaps across replicates
olap.rate.control <- dba.overlap(IMR32_data,IMR32_data$masks$control,mode=DBA_OLAP_RATE) 
print(olap.rate.control)
olap.rate.PB <- dba.overlap(IMR32_data,IMR32_data$masks$PB5d,mode=DBA_OLAP_RATE) 
print(olap.rate.PB)

# check peak overlaps across all samples
olap.rate.all <- dba.overlap(IMR32_data,mode=DBA_OLAP_RATE) 
print(olap.rate.all) 

# get consensus peak sets for each condition
#keep if in at least 3 of the 5 control replicates
control_consensus_1<-dba.peakset(IMR32_data, bRetrieve=T, DataType=DBA_DATA_FRAME,IMR32_data$masks$control, minOverlap=3) 
control_consensus_1x<-control_consensus_1[,1:3]
nrow(control_consensus_1x) 
write.table(control_consensus_1x,paste0(datadir,'PB_H3K27ac_ChIPseq/output_data/IMR32_control_ConsensusNoSummits_vInput_broad_in3of5.bed'), col.names = F, row.names = F, quote = F, sep = "\t")

PB_consensus_1<-dba.peakset(IMR32_data, bRetrieve=T, DataType=DBA_DATA_FRAME,IMR32_data$masks$PB, minOverlap=3) 
PB_consensus_1x<-PB_consensus_1[,1:3]
nrow(PB_consensus_1x) 
write.table(PB_consensus_1x,paste0(datadir,'PB_H3K27ac_ChIPseq/output_data/IMR32_PB_ConsensusNoSummits_vInput_broad_in3of5.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


# get full consensus peakset across conditions
IMR32_consensus <- dba.peakset(IMR32_data, consensus = DBA_CONDITION, minOverlap=3)
ExptConsensus <- dba(IMR32_consensus, mask=IMR32_consensus$masks$Consensus,minOverlap=1)
ExptConsensus
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE,DataType=DBA_DATA_FRAME,minOverlap=1)
ConsensusPeaksx<-ConsensusPeaks[,1:3]
nrow(ConsensusPeaksx) 
write.table(ConsensusPeaksx,paste0(datadir,'PB_H3K27ac_ChIPseq/output_data/IMR32_ConsensusNoSummits_vInput_broad_bothin3of5.bed'), col.names = F, row.names = F, quote = F, sep = "\t")


# get peakset counts ------------------------------------------------------
register(MulticoreParam(4))
ConsensusPeaks <- dba.peakset(ExptConsensus, bRetrieve=TRUE)
IMR32_data_count <- dba.count(IMR32_data,peaks=ConsensusPeaks,bParallel=TRUE,summits=FALSE) 


# plot PCA ----------------------------------------------------------------
tiff(file='PB_H3K27ac_ChIPseq/output_figures/IMR32_PCA_normalised_ConsensusNoSummits_vInput_broad_nolabel.tiff',width=1400,height=1300,res=300)
dba.plotPCA(IMR32_data_count_v2,DBA_CONDITION,vColors=c('grey48','red3'),score=DBA_SCORE_NORMALIZED)
dev.off()


# get FRiP details; FRiP = fraction of reads in the peaks -----------------
info <- dba.show(IMR32_data_count_v2)
write.table(data.frame(info,stringsAsFactors=FALSE),'PB_H3K27ac_ChIPseq/output_figures/info_IMR32_vInput_broad.txt',col.names=TRUE,row.names=FALSE,quote=FALSE)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
print(libsizes)


# normalisation -----------------------------------------------------------
IMR32_data_normalize<-dba.normalize(IMR32_data_count_v2,method=DBA_DESEQ2) 


# conduct differential analysis -------------------------------------------

# sort contrast
IMR32_data_contrast <- dba.contrast(IMR32_data_normalize,design="~Replicate+Condition") 

# run analysis
IMR32_data_analyze <- dba.analyze(IMR32_data_contrast,method = DBA_DESEQ2,bGreylist=FALSE) 
dba.show(IMR32_data_analyze, bContrasts=TRUE)

# check differential sites
IMR32_data_analyze_DiffSites <- dba.report(IMR32_data_analyze,method = DBA_DESEQ2) 
print(sum(IMR32_data_analyze_DiffSites$Fold>0 & IMR32_data_analyze_DiffSites$FDR<0.05))
print(sum(IMR32_data_analyze_DiffSites$Fold<0 & IMR32_data_analyze_DiffSites$FDR<0.05))



# write outputs -----------------------------------------------------------

# DiffBind output table
report <- dba.report(IMR32_data_analyze, contrast = 1, th = 1, bFlip = F)
write.table(report, paste0("PB_H3K27ac_ChIPseq/output_data/IMR32_5dPBvControl_ConsensusNoSummits_vInput.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

# plot data as volcano plot
png('PB_H3K27ac_ChIPseq/output_figures/IMR32_volcano_afteranalysis_ConsensusNoSummits_vInput_broad.png',width=500,height=400)
dba.plotVolcano(IMR32_data_analyze) 
dev.off()

report_df<-data.frame(report)
nrow(report_df)
nrow(subset(report_df,FDR<0.05 & Fold>0.5))
nrow(subset(report_df,FDR<0.05 & Fold<(-0.5)))
nrow(subset(report_df,FDR>0.5))



# get sub-group files for downstream plotting ---------------

# for input into ChIPseeker
write.table(subset(report_df,FDR<0.05 & Fold>0.5), "PB_H3K27ac_ChIPseq/output_data/IMR32_5dPBvControl_UP_FDR005Fold05.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(subset(report_df,FDR<0.05 & Fold<(-0.5)), "PB_H3K27ac_ChIPseq/output_data/IMR32_5dPBvControl_DOWN_FDR005Fold05.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(subset(report_df,FDR>0.5), "PB_H3K27ac_ChIPseq/output_data/IMR32_5dPBvControl_ns_FDR05.txt", sep = "\t", col.names = T, row.names = F, quote = F)


# for input into plotHeatmap
IMR32_broad<-report_df
IMR32_broad$name<-apply(IMR32_broad[,1:3],1,function(x) paste(x[1:3],collapse=':'))
IMR32_broad$name<-gsub(' ','',IMR32_broad$name)

#order outputs by FDR
IMR32_broad<-IMR32_broad[order(IMR32_broad$FDR,decreasing=FALSE),]
IMR32_broad_sigup<-subset(IMR32_broad,Fold>0.5 & FDR<0.05)
nrow(IMR32_broad_sigup)
IMR32_broad_sigdown<-subset(IMR32_broad,Fold<(-0.5) & FDR<0.05)
nrow(IMR32_broad_sigdown)
IMR32_broad_ns<-subset(IMR32_broad,FDR>0.5)
nrow(IMR32_broad_ns)

write.table(IMR32_broad_sigup[,c(1,2,3,4,5,11)],'PB_H3K27ac_ChIPseq/output_data/IMR32_broadpeak_up05.bed',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(IMR32_broad_sigdown[,c(1,2,3,4,5,11)],'PB_H3K27ac_ChIPseq/output_data/IMR32_broadpeak_down05.bed',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(IMR32_broad_ns[,c(1,2,3,4,5,11)],'PB_H3K27ac_ChIPseq/output_data/IMR32_broadpeak_ns05.bed',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)



