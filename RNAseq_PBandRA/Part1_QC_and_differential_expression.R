#this script conducts data filtering and differential expression analysis with DESeq2 for the palbociclib and retinoic acid RNA-seq data


# load libraries ----------------------------------------------------------
library(DESeq2)
library(limma)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# read in the data and sample information --------------------------------------------------------

# raw counts
RNAseq_counts<-read.csv('PBandRA_RNAseq/data/BE2C_PBandRA_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
names(RNAseq_counts)
dim(RNAseq_counts)

# sample details
sample_info<-read.csv('PBandRA_RNAseq/data/sample_information.csv',stringsAsFactors=FALSE,header=TRUE)

# sort format of sample details
sample_info$sample<-apply(sample_info,1,function(x) paste(x[(length(x)-1)],x[length(x)],sep='_'))
sample_info<-sample_info[,c('sample','condition','replicate')]
sample_info$condition<-factor(sample_info$condition,levels=c('DMSO','PB','RA','PBandRA'))
sample_info<-sample_info[order(sample_info$condition),]
sample_info



# generate correlation matrix and PCA plot ----------------------------------------------------------------

counts_data<-RNAseq_counts[rowSums(RNAseq_counts[,2:17])>10,]
names(counts_data)
dim(counts_data)
rownames(counts_data)<-counts_data$gene_name
counts_data<-counts_data[,-1]

dds <- DESeqDataSetFromMatrix(countData=counts_data,colData=sample_info,design=~replicate+condition)

# filtering
#dds <- dds[ rowSums(counts(dds)) > 10, ] 

# normalisation and preprocessing
dds <- DESeq(dds)

# VST normalisation
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd))) 

# Correlation plot of all samples
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Replicate, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Condition, vsd$Replicate, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

tiff('PBandRA_RNAseq/output_figures/vst_corrPlot_all.tiff', height = 2000, width = 2200,res=300)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
dev.off()


# plot PCAs for different condition combinations
generate_PCA <- function(vsd,groups,colourset,setname) {
  vsd.sub <- vsd[ , vsd$condition %in% groups]
  pcaData <- DESeq2::plotPCA(vsd.sub, intgroup = c( "replicate", "condition"), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
  percentVar <- round(100 * attr(pcaData, "percentVar")) 
  
  f1<-ggplot(pcaData, aes(x = PC1, y = PC2, fill = factor(Condition),col = factor(Condition), shape = factor(Replicate))) + geom_point(size =3) + 
    scale_shape_manual(values=c(21,22,23,24)) + scale_fill_manual(values=colourset) + scale_colour_manual(values=colourset) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text=element_text(size=12))
  ggsave(plot=f1,file=paste0('PBandRA_RNAseq/output_figures/',setname,"_PCA.tiff"),width=3.75,height=2.5)
  
}

setnames<-c('DMSOvPB','DMSOvRA','DMSOvPBRA','all')
groupsets<-list(c('DMSO','PB'),c('DMSO','RA'),c('DMSO','PBandRA'),c('DMSO','RA','PB','PBandRA'))
coloursets<-list(c('grey58','red3'),c('grey58','dodgerblue3'),c('grey58','purple3'),c('grey58','red3','dodgerblue3','purple3'))

for(i in 1:length(groupsets)){
  generate_PCA(vsd,groupsets[[i]],coloursets[[i]],setnames[i])
}



# differential expression with DESeq2 -------------------------------------

## for all comparisons to control: PB v Control, RA v Control, PB+RA v Control ##

# Select data
counts_matrix<-subset(RNAseq_counts,rowSums(RNAseq_counts[,2:ncol(RNAseq_counts)],na.rm=TRUE)>10)
dim(counts_matrix)
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(ccounts_matrix[,c(-1)])
head(counts_matrix)

# do the analysis
DESeq2data<-DESeqDataSetFromMatrix(countData = countsMatrix,colData = sample_info,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 

#PB v Control
resLFC <- lfcShrink(DESeq2output, coef='condition_PB_vs_DMSO', type="apeglm") 
#RA v Control
resLFC2 <- lfcShrink(DESeq2output, coef='condition_RA_vs_DMSO', type="apeglm") 
#PBandRA v Control
resLFC3 <- lfcShrink(DESeq2output, coef='condition_PBandRA_vs_DMSO', type="apeglm") 

rm(DESeq2data,DESeq2output)

# quick look at results
png("PBandRA_RNAseq/output_figures/MAplot_PBvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PBandRA_RNAseq/output_figures/MAplot_RAvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC2, ylim=c(-4,4))
dev.off()

png("PBandRA_RNAseq/output_figures/MAplot_PBandRAvControl.png",width=450,height=400)
DESeq2::plotMA(resLFC3, ylim=c(-4,4)) 
dev.off()

# results table formatting
PBvControl_results<-as.data.frame(resOrdered)
nrow(PBvControl_results)
PBvControl_results$gene<-rownames(PBvControl_results)
PBvControl_results<-PBvControl_results[,c(6,1:5)]

RAvControl_results<-as.data.frame(resOrdered2)
nrow(RAvControl_results)
RAvControl_results$gene<-rownames(RAvControl_results)
RAvControl_results<-RAvControl_results[,c(6,1:5)]

PBandRAvControl_results<-as.data.frame(resOrdered3)
nrow(PBandRAvControl_results)
PBandRAvControl_results$gene<-rownames(PBandRAvControl_results)
PBandRAvControl_results<-PBandRAvControl_results[,c(6,1:5)]

# write results
write.csv(PBvControl_results, file="PBandRA_RNAseq/output_data/PBvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(RAvControl_results, file="PBandRA_RNAseq/output_data/RAvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(PBandRAvControl_results, file="PBandRA_RNAseq/output_data/PBandRAvControl_DESeq2output.csv",row.names=FALSE,quote=FALSE)

rm(resLFC,resLFC2,resLFC3)
rm(resOrdered,resOrdered2,resOrdered3)




## for comparison: PB+RA v PB ##

# get subset of sample annotation and counts matrix
counts_matrix<-subset(RNAseq_counts,rowSums(RNAseq_counts[,2:ncol(RNAseq_counts)],na.rm=TRUE)>10)
dim(counts_matrix)
names(counts_matrix)
counts_matrix<-counts_matrix[,c(1,6:17)]
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

sample_info<-read.csv('PBandRA_RNAseq/data/sample_information.csv',header=TRUE,stringsAsFactors = FALSE)
sample_info
sample_info$sample<-apply(sample_info,1,function(x) paste(x[(length(x)-1)],x[length(x)],sep='_'))
sample_info<-sample_info[,c('sample','condition','replicate')]
sample_info
sample_info$condition<-factor(sample_info$condition,levels=c('DMSO','PB','RA','PBandRA'))
sample_info<-sample_info[order(sample_info$condition),]
sample_info
sample_info<-sample_info[5:16,]
sample_info

# do the analysis
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = sample_info,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 
#PBandRA v PB
resLFC <- lfcShrink(DESeq2output, coef='condition_PBandRA_vs_PB', type="apeglm") 
rm(DESeq2data,DESeq2output)

# sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

# quick look at results
png("PBandRA_RNAseq/output_figures/MAplot_PBandRAvPB.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

# results table
PBandRAvPB_results<-as.data.frame(resOrdered)
nrow(PBandRAvPB_results)
PBandRAvPB_results$gene<-rownames(PBandRAvPB_results)
PBandRAvPB_results<-PBandRAvPB_results[,c(6,1:5)]

# write results
write.csv(PBandRAvPB_results, file="PBandRA_RNAseq/output_data/PBandRAvPB_DESeq2output.csv",row.names=FALSE,quote=FALSE)

rm(resLFC)
rm(resOrdered)




## for comparisons: PB+RA v RA, PB v RA ##

# get subset of sample annotation and counts matrix
counts_matrix<-subset(RNAseq_counts,rowSums(RNAseq_counts[,2:ncol(RNAseq_counts)],na.rm=TRUE)>10)
dim(counts_matrix)
names(counts_matrix)
counts_matrix<-counts_matrix[,c(1,6:17)]
names(counts_matrix)
rownames(counts_matrix)<-counts_matrix[,1]
counts_matrix<-as.matrix(counts_matrix[,c(-1)])
head(counts_matrix)

sample_info<-read.csv('PBandRA_RNAseq/data/sample_information.csv',header=TRUE,stringsAsFactors = FALSE)
sample_info
sample_info$sample<-apply(sample_info,1,function(x) paste(x[(length(x)-1)],x[length(x)],sep='_'))
sample_info<-sample_info[,c('sample','condition','replicate')]
sample_info
sample_info$condition<-factor(sample_info$condition,levels=c('DMSO','PB','RA','PBandRA'))
sample_info<-sample_info[order(sample_info$condition),]
sample_info
sample_info<-sample_info[5:16,]
sample_info

sample_info$condition<-factor(sample_info$condition,levels=c('RA','PB','PBandRA'))
DESeq2data<-DESeqDataSetFromMatrix(countData = counts_matrix,colData = sample_info,design= ~ replicate + condition)
DESeq2output <- DESeq(DESeq2data)
resultsNames(DESeq2output) 

#PBandRA v RA
resLFC <- lfcShrink(DESeq2output, coef='condition_PBandRA_vs_RA', type="apeglm") 
#PB v RA
resLFC2 <- lfcShrink(DESeq2output, coef='condition_PB_vs_RA', type="apeglm") 
rm(DESeq2data,DESeq2output)

# sort results
resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resLFC)
sum(resLFC$padj < 0.05, na.rm=TRUE) 

resOrdered2 <- resLFC2[order(resLFC2$pvalue),]
summary(resLFC2)
sum(resLFC2$padj < 0.05, na.rm=TRUE) 


# quick look at results
png("PBandRA_RNAseq/output_figures/MAplot_PBandRAvRA.png",width=450,height=400)
DESeq2::plotMA(resLFC, ylim=c(-4,4))
dev.off()

png("PBandRA_RNAseq/output_figures/MAplot_PBvRA.png",width=450,height=400)
DESeq2::plotMA(resLFC2, ylim=c(-4,4))
dev.off()

# results table
PBandRAvRA_results<-as.data.frame(resOrdered)
nrow(PBandRAvRA_results)
PBandRAvRA_results$gene<-rownames(PBandRAvRA_results)
PBandRAvRA_results<-PBandRAvRA_results[,c(6,1:5)]

PBvRA_results<-as.data.frame(resOrdered2)
nrow(PBvRA_results)
PBvRA_results$gene<-rownames(PBvRA_results)
PBvRA_results<-PBvRA_results[,c(6,1:5)]

# write results
write.csv(PBandRAvRA_results, file="PBandRA_RNAseq/output_data/PBandRAvRA_DESeq2output.csv",row.names=FALSE,quote=FALSE)
write.csv(PBvRA_results, file="PBandRA_RNAseq/output_data/differential_expression/PBvRA_DESeq2output.csv",row.names=FALSE,quote=FALSE)

rm(resLFC,resLFC2)
rm(resOrdered,resOrdered2)











