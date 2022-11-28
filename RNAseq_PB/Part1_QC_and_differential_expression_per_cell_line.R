#this script conducts data filtering and differential expression analysis with DESeq2 for the palbociclib RNA-seq data in three neuroblastoma cell lines


# load libraries ----------------------------------------------------------
library(DESeq2)
library(limma)
library(apeglm)
library(ggplot2)



########################################################################################
################################## SK-N-BE(2)C data ####################################
########################################################################################

# read in the data and sample information --------------------------------------------------------

# raw counts
BE2C_RNAseq_counts<-read.csv('PB_RNAseq/data/BE2C_PB_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
names(BE2C_RNAseq_counts)
dim(BE2C_RNAseq_counts)

# sample details
BE2C_sampleinfo<-read.csv('PB_RNAseq/data/sampleinfo_BE2C_PB_RNAseq.csv',stringsAsFactors = FALSE,header=TRUE)
BE2C_sampleinfo
BE2C_sampleinfo$condition<-factor(BE2C_sampleinfo$condition,levels=c('control','24hPB','7dPB'))
BE2C_sampleinfo$replicate<-factor(BE2C_sampleinfo$replicate,levels=c('R1','R2','R3','R4','R5'))


# PCA plot ----------------------------------------------------------------

# BE2C
count_df<-BE2C_RNAseq_counts
names(count_df)
dim(count_df)
rownames(count_df)<-BE2C_RNAseq_counts$gene
count_df<-count_df[rowSums(count_df[,2:ncol(count_df)])>10,]
count_df<-count_df[,-1]
DESeq2data <- DESeqDataSetFromMatrix(countData=as.matrix(count_df),colData=BE2C_sampleinfo,design=~replicate+condition)
DESeq2output <- DESeq(DESeq2data)
vsd <- vst(DESeq2output, blind=FALSE)
pca_data<-plotPCA(vsd, intgroup=c("condition", "replicate"),returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar")) 
f1<-ggplot(pca_data, aes(x = PC1, y = PC2, fill = factor(condition),col = factor(condition), shape = factor(replicate))) + geom_point(size =3) + 
  scale_shape_manual(values=c(21,22,23,24,25)) + scale_fill_manual(values=c('grey58','goldenrod2','red3')) + scale_colour_manual(values=c('grey58','goldenrod2','red3')) + 
  xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) + theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave(plot=f1,file=paste0("PB_RNAseq/output_figures/BE2C_PCA.png"),width=4,height=2.6)

rm(count_df,DESeq2data,DESeq2output,vsd,pca_data,percent_var,f1)


# differential expression with DESeq2 -------------------------------------

# select data
counts_matrix<-subset(BE2C_RNAseq_counts,rowSums(BE2C_RNAseq_counts[,2:ncol(BE2C_RNAseq_counts)],na.rm=TRUE)>30)
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

# do the analysis
#24hPB v control
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = BE2C_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC <- lfcShrink(DESeq2output, coef='condition_24hPB_vs_control', type="apeglm") 

rm(DESeq2data,DESeq2output)

#7dPB v control
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = BE2C_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC2 <- lfcShrink(DESeq2output, coef='condition_7dPB_vs_control', type="apeglm") 

rm(DESeq2data,DESeq2output)

#7dPB v 24hPB
BE2C_sampleinfo<-BE2C_sampleinfo[grepl('PB',BE2C_sampleinfo$condition)==TRUE,]
BE2C_sampleinfo

counts_matrix<-subset(BE2C_RNAseq_counts,rowSums(BE2C_RNAseq_counts[,2:ncol(BE2C_RNAseq_counts)],na.rm=TRUE)>30)[,grepl('gene|PB',names(BE2C_RNAseq_counts))==TRUE]
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = BE2C_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC3 <- lfcShrink(DESeq2output, coef='condition_7dPB_vs_24hPB', type="apeglm") 

rm(DESeq2data,DESeq2output)

# sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

resOrdered2 <- resLFC2[order(resLFC2$pvalue),]
summary(resLFC2)
sum(resLFC2$padj < 0.05, na.rm=TRUE) 

resOrdered3 <- resLFC3[order(resLFC3$pvalue),]
summary(resLFC3)
sum(resLFC3$padj < 0.05, na.rm=TRUE) 

# quick look at results
png("PB_RNAseq/output_figures/MAplot_24hPBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PB_RNAseq/output_figures/MAplot_7dPBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PB_RNAseq/output_figures/MAplot_7dPBv24hPB.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()


# results table
BE2C_24hPBvControl_results<-as.data.frame(resOrdered)
nrow(BE2C_24hPBvControl_results)
BE2C_24hPBvControl_results$gene<-rownames(BE2C_24hPBvControl_results)
BE2C_24hPBvControl_results<-BE2C_24hPBvControl_results[,c(6,1:5)]

BE2C_7dPBvControl_results<-as.data.frame(resOrdered2)
nrow(BE2C_7dPBvControl_results)
BE2C_7dPBvControl_results$gene<-rownames(BE2C_7dPBvControl_results)
BE2C_7dPBvControl_results<-BE2C_7dPBvControl_results[,c(6,1:5)]

BE2C_7dPBv24hPB_results<-as.data.frame(resOrdered3)
nrow(BE2C_7dPBv24hPB_results)
BE2C_7dPBv24hPB_results$gene<-rownames(BE2C_7dPBv24hPB_results)
BE2C_7dPBv24hPB_results<-BE2C_7dPBv24hPB_results[,c(6,1:5)]

# write results
write.csv(BE2C_24hPBvControl_results, file="PB_RNAseq/output_data/BE2C_24hPBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(BE2C_7dPBvControl_results, file="PB_RNAseq/output_data/BE2C_7dPBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(BE2C_7dPBv24hPB_results, file="PB_RNAseq/output_data/BE2C_7dPBv24hPB_DESeq2output.csv",row.names=FALSE,quote=FALSE)

rm(resLFC,resLFC2,resLFC3)
rm(resOrdered,resOrdered2,resOrdered3)




#######################################################################################
#################################### IMR-32 data ######################################
#######################################################################################

# read in the data and sample information --------------------------------------------------------

# raw counts
IMR32_RNAseq_counts<-read.csv('PB_RNAseq/data/IMR32_PB_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
names(IMR32_RNAseq_counts)
dim(IMR32_RNAseq_counts)

# sample details
IMR32_sampleinfo<-read.csv('PB_RNAseq/data/sampleinfo_IMR32_PB_RNAseq.csv',stringsAsFactors = FALSE,header=TRUE)
IMR32_sampleinfo
IMR32_sampleinfo$condition<-factor(IMR32_sampleinfo$condition,levels=c('control','24hPB','5dPB'))
IMR32_sampleinfo$replicate<-factor(IMR32_sampleinfo$replicate,levels=c('R1','R2','R3','R4','R5'))


# PCA plot ----------------------------------------------------------------

# IMR32
count_df<-IMR32_RNAseq_counts
names(count_df)
dim(count_df)
rownames(count_df)<-IMR32_RNAseq_counts$gene
count_df<-count_df[rowSums(count_df[,2:ncol(count_df)])>10,]
count_df<-count_df[,-1]
DESeq2data <- DESeqDataSetFromMatrix(countData=as.matrix(count_df),colData=IMR32_sampleinfo,design=~replicate+condition)
DESeq2output <- DESeq(DESeq2data)
vsd <- vst(DESeq2output, blind=FALSE)
pca_data<-plotPCA(vsd, intgroup=c("condition", "replicate"),returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar")) 
f1<-ggplot(pca_data, aes(x = PC1, y = PC2, fill = factor(condition),col = factor(condition), shape = factor(replicate))) + geom_point(size =3) + 
  scale_shape_manual(values=c(21,22,23,24,25)) + scale_fill_manual(values=c('grey58','goldenrod2','red3')) + scale_colour_manual(values=c('grey58','goldenrod2','red3')) + 
  xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) + theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave(plot=f1,file=paste0("PB_RNAseq/output_figures/IMR32_PCA.png"),width=4,height=2.6)

rm(count_df,DESeq2data,DESeq2output,vsd,pca_data,percent_var,f1)


# differential expression with DESeq2 -------------------------------------

# select data
counts_matrix<-subset(IMR32_RNAseq_counts,rowSums(IMR32_RNAseq_counts[,2:ncol(IMR32_RNAseq_counts)],na.rm=TRUE)>30)
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

# do the analysis
#24hPB v control
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = IMR32_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC <- lfcShrink(DESeq2output, coef='condition_24hPB_vs_control', type="apeglm") 

rm(DESeq2data,DESeq2output)

#5dPB v control
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = IMR32_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC2 <- lfcShrink(DESeq2output, coef='condition_5dPB_vs_control', type="apeglm") 

rm(DESeq2data,DESeq2output)

#5dPB v 24hPB
IMR32_sampleinfo<-IMR32_sampleinfo[grepl('PB',IMR32_sampleinfo$condition)==TRUE,]
IMR32_sampleinfo

counts_matrix<-subset(IMR32_RNAseq_counts,rowSums(IMR32_RNAseq_counts[,2:ncol(IMR32_RNAseq_counts)],na.rm=TRUE)>30)[,grepl('gene|PB',names(IMR32_RNAseq_counts))==TRUE]
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = IMR32_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC3 <- lfcShrink(DESeq2output, coef='condition_5dPB_vs_24hPB', type="apeglm") 

rm(DESeq2data,DESeq2output)

# sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

resOrdered2 <- resLFC2[order(resLFC2$pvalue),]
summary(resLFC2)
sum(resLFC2$padj < 0.05, na.rm=TRUE) 

resOrdered3 <- resLFC3[order(resLFC3$pvalue),]
summary(resLFC3)
sum(resLFC3$padj < 0.05, na.rm=TRUE) 

#quick look at results
png("PB_RNAseq/output_figures/MAplot_24hPBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PB_RNAseq/output_figures/MAplot_5dPBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PB_RNAseq/output_figures/MAplot_5dPBv24hPB.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()


# results table
IMR32_24hPBvControl_results<-as.data.frame(resOrdered)
nrow(IMR32_24hPBvControl_results)
IMR32_24hPBvControl_results$gene<-rownames(IMR32_24hPBvControl_results)
IMR32_24hPBvControl_results<-IMR32_24hPBvControl_results[,c(6,1:5)]

IMR32_5dPBvControl_results<-as.data.frame(resOrdered2)
nrow(IMR32_5dPBvControl_results)
IMR32_5dPBvControl_results$gene<-rownames(IMR32_5dPBvControl_results)
IMR32_5dPBvControl_results<-IMR32_5dPBvControl_results[,c(6,1:5)]

IMR32_5dPBv24hPB_results<-as.data.frame(resOrdered3)
nrow(IMR32_5dPBv24hPB_results)
IMR32_5dPBv24hPB_results$gene<-rownames(IMR32_5dPBv24hPB_results)
IMR32_5dPBv24hPB_results<-IMR32_5dPBv24hPB_results[,c(6,1:5)]

# write results
write.csv(IMR32_24hPBvControl_results, file="PB_RNAseq/output_data/IMR32_24hPBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(IMR32_5dPBvControl_results, file="PB_RNAseq/output_data/IMR32_5dPBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(IMR32_5dPBv24hPB_results, file="PB_RNAseq/output_data/IMR32_5dPBv24hPB_DESeq2output.csv",row.names=FALSE,quote=FALSE)

rm(resLFC,resLFC2,resLFC3)
rm(resOrdered,resOrdered2,resOrdered3)


########################################################################################
#################################### SH-SY5Y data ######################################
########################################################################################

# read in the data and sample information --------------------------------------------------------

# raw counts
SHSY5Y_RNAseq_counts<-read.csv('PB_RNAseq/data/SHSY5Y_PB_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
names(SHSY5Y_RNAseq_counts)
dim(SHSY5Y_RNAseq_counts)

# sample details
SHSY5Y_sampleinfo<-read.csv('PB_RNAseq/data/sampleinfo_SHSY5Y_PB_RNAseq.csv',stringsAsFactors = FALSE,header=TRUE)
SHSY5Y_sampleinfo
SHSY5Y_sampleinfo$condition<-factor(SHSY5Y_sampleinfo$condition,levels=c('control','24hPB','5dPB'))
SHSY5Y_sampleinfo$replicate<-factor(SHSY5Y_sampleinfo$replicate,levels=c('R1','R2','R3','R4','R5'))


# PCA plot ----------------------------------------------------------------

# SHSY5Y
count_df<-SHSY5Y_RNAseq_counts
names(count_df)
dim(count_df)
rownames(count_df)<-SHSY5Y_RNAseq_counts$gene
count_df<-count_df[rowSums(count_df[,2:ncol(count_df)])>10,]
count_df<-count_df[,-1]
DESeq2data <- DESeqDataSetFromMatrix(countData=as.matrix(count_df),colData=SHSY5Y_sampleinfo,design=~replicate+condition)
DESeq2output <- DESeq(DESeq2data)
vsd <- vst(DESeq2output, blind=FALSE)
pca_data<-plotPCA(vsd, intgroup=c("condition", "replicate"),returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar")) 
f1<-ggplot(pca_data, aes(x = PC1, y = PC2, fill = factor(condition),col = factor(condition), shape = factor(replicate))) + geom_point(size =3) + 
  scale_shape_manual(values=c(21,22,23,24,25)) + scale_fill_manual(values=c('grey58','goldenrod2','red3')) + scale_colour_manual(values=c('grey58','goldenrod2','red3')) + 
  xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) + theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave(plot=f1,file=paste0("PB_RNAseq/output_figures/SHSY5Y_PCA.png"),width=4,height=2.6)

rm(count_df,DESeq2data,DESeq2output,vsd,pca_data,percent_var,f1)


# differential expression with DESeq2 -------------------------------------

# select data
counts_matrix<-subset(SHSY5Y_RNAseq_counts,rowSums(SHSY5Y_RNAseq_counts[,2:ncol(SHSY5Y_RNAseq_counts)],na.rm=TRUE)>30)
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

# do the analysis
#24hPB v control
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = SHSY5Y_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC <- lfcShrink(DESeq2output, coef='condition_24hPB_vs_control', type="apeglm") 

rm(DESeq2data,DESeq2output)

#5dPB v control
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = SHSY5Y_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC2 <- lfcShrink(DESeq2output, coef='condition_5dPB_vs_control', type="apeglm") 

rm(DESeq2data,DESeq2output)

#5dPB v 24hPB
SHSY5Y_sampleinfo<-SHSY5Y_sampleinfo[grepl('PB',SHSY5Y_sampleinfo$condition)==TRUE,]
SHSY5Y_sampleinfo

counts_matrix<-subset(SHSY5Y_RNAseq_counts,rowSums(SHSY5Y_RNAseq_counts[,2:ncol(SHSY5Y_RNAseq_counts)],na.rm=TRUE)>30)[,grepl('gene|PB',names(SHSY5Y_RNAseq_counts))==TRUE]
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = SHSY5Y_sampleinfo,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
resLFC3 <- lfcShrink(DESeq2output, coef='condition_5dPB_vs_24hPB', type="apeglm") 

rm(DESeq2data,DESeq2output)

# sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

resOrdered2 <- resLFC2[order(resLFC2$pvalue),]
summary(resLFC2)
sum(resLFC2$padj < 0.05, na.rm=TRUE) 

resOrdered3 <- resLFC3[order(resLFC3$pvalue),]
summary(resLFC3)
sum(resLFC3$padj < 0.05, na.rm=TRUE) 

# quick look at results
png("PB_RNAseq/output_figures/MAplot_24hPBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PB_RNAseq/output_figures/MAplot_5dPBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PB_RNAseq/output_figures/MAplot_5dPBv24hPB.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()


# results table
SHSY5Y_24hPBvControl_results<-as.data.frame(resOrdered)
nrow(SHSY5Y_24hPBvControl_results)
SHSY5Y_24hPBvControl_results$gene<-rownames(SHSY5Y_24hPBvControl_results)
SHSY5Y_24hPBvControl_results<-SHSY5Y_24hPBvControl_results[,c(6,1:5)]

SHSY5Y_5dPBvControl_results<-as.data.frame(resOrdered2)
nrow(SHSY5Y_5dPBvControl_results)
SHSY5Y_5dPBvControl_results$gene<-rownames(SHSY5Y_5dPBvControl_results)
SHSY5Y_5dPBvControl_results<-SHSY5Y_5dPBvControl_results[,c(6,1:5)]

SHSY5Y_5dPBv24hPB_results<-as.data.frame(resOrdered3)
nrow(SHSY5Y_5dPBv24hPB_results)
SHSY5Y_5dPBv24hPB_results$gene<-rownames(SHSY5Y_5dPBv24hPB_results)
SHSY5Y_5dPBv24hPB_results<-SHSY5Y_5dPBv24hPB_results[,c(6,1:5)]

# write results
write.csv(SHSY5Y_24hPBvControl_results, file="PB_RNAseq/output_data/SHSY5Y_24hPBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(SHSY5Y_5dPBvControl_results, file="PB_RNAseq/output_data/SHSY5Y_5dPBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(SHSY5Y_5dPBv24hPB_results, file="PB_RNAseq/output_data/SHSY5Y_5dPBv24hPB_DESeq2output.csv",row.names=FALSE,quote=FALSE)

rm(resLFC,resLFC2,resLFC3)
rm(resOrdered,resOrdered2,resOrdered3)








