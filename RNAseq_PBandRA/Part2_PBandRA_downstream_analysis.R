#this script conducts downstream investigations of the palbociclib and retinoic acid RNA-seq data


# load libraries ----------------------------------------------------------

library(FSA)
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(factoextra)
library(scales)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

# select significantly changed gene set -----------------------------------

#combine differential expression tables from the different comparisons conducted
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


data_path='PBandRA_RNAseq/output_data/'
filenames <- dir(data_path, pattern = "*.csv")
length(filenames)

allinputs<-lapply(as.list(filenames),function(x) read.delim(paste0(data_path,x),sep=',',stringsAsFactors = FALSE,header=TRUE))
length(allinputs)

allDE_collated<-collate_results(allinputs, filenames)
names(allDE_collated)
dim(allDE_collated)

# select all genes that have a significant change in at least one condition with a log2FC > 0.5 or < -0.5
alldiffgeneset<-subset(allDE_collated,PBandRAvControl_padj<0.05 | PBandRAvPB_padj<0.05 | PBandRAvRA_padj<0.05 | PBvControl_padj<0.05 | RAvControl_padj<0.05)
nrow(alldiffgeneset) 
alldiffgeneset<-subset(alldiffgeneset, PBandRAvControl_log2FoldChange>0.5 | PBandRAvControl_log2FoldChange<(-0.5) | PBandRAvPB_log2FoldChange>0.5 | PBandRAvPB_log2FoldChange<(-0.5) | PBandRAvRA_log2FoldChange>0.5 | PBandRAvRA_log2FoldChange<(-0.5) | PBvControl_log2FoldChange>0.5 | PBvControl_log2FoldChange<(-0.5) | RAvControl_log2FoldChange>0.5 | RAvControl_log2FoldChange<(-0.5))
nrow(alldiffgeneset) 


# get CPM normalised counts -----------------------------------------------

#raw counts
RNAseq_counts<-read.csv('PBandRA_RNAseq/data/BE2C_PBandRA_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
names(RNAseq_counts)
dim(RNAseq_counts)

#convert to CPM (counts per million)
RNAseq_CPM<-RNAseq_counts
RNAseq_CPM[,2:17]<-apply(RNAseq_CPM[,2:17],2,function(x) (x/sum(x))*1000000)
lapply(RNAseq_CPM[,2:17],function(x) sum(x))

#write output table
write.csv(RNAseq_CPM,file='PBandRA_RNAseq/output_data/BE2C_PBandRA_RNAseq_CPM.csv',quote=FALSE,row.names=FALSE)

 
# conduct kmeans clustering -----------------------------------------------

#get normalised counts for all genes that have are significantly differential in at least one comparison
diffgene_CPM<-subset(RNAseq_CPM,(gene_name %in% alldiffgeneset$gene)==TRUE)
nrow(diffgene_CPM)

#row scaling of CPM normalised counts
diffgenes_prescale<-diffgene_CPM
for(i in 1:nrow(diffgenes_prescale)){
  diffgenes_prescale[i,2:ncol(diffgenes_prescale)]<-scale(as.numeric(diffgenes_prescale[i,2:ncol(diffgenes_prescale)]))
}

kdf<-diffgenes_prescale
rownames(kdf)<-diffgenes_prescale$gene_name
nrow(kdf)
rownames(kdf)<-kdf$gene_name
kdf<-kdf[,-1]

#determine optimal number of clusters
#NOTE: with the selection of significantly differentially expressed either up or down regulated genes here, two clusters will dominate first
png('PBandRA_RNAseq/output_figures/kmeans_clusters_elbowmethod.png',width=400,height=250)
fviz_nbclust(kdf, kmeans, method='wss')
dev.off()

png('PBandRA_RNAseq/output_figures/kmeans_clusters_elbowmethod_zoom.png',width=400,height=250)
fviz_nbclust(kdf, kmeans, method='wss')+coord_trans(ylim=c(25000,75000))
dev.off()

#visual check of the value of k choice
res_k<-hkmeans(kdf,5,hc.method='centroid')
png('PBandRA_RNAseq/output_figures/kmeans_clusters_5.png',width=350,height=300)
fviz_cluster(res_k, ellipse=FALSE,labelsize=6,ggtheme = theme_classic(),geom=c('point'),pointsize=1,stand=TRUE,shape='circle')+scale_colour_manual(values = c("darkorange", "green4","mediumvioletred",'royalblue3','darkgoldenrod2','red3','grey78','grey28'))
dev.off()

#conduct kmeans clustering with k=5
res_k<-hkmeans(kdf,5,hc.method='complete')
kdf$hk_clusters<-res.k$cluster
kdf$hk_clusters<-as.character(kdf$hk_clusters)
kdf<-kdf[order(kdf$hk_clusters),]

#extract number of genes in each cluster
for(i in 1:5){
  print(nrow(subset(kdf,hk_clusters==i)))
}      
    
       
# plot heatmap of clusters ------------------------------------------------

h1<-Heatmap(as.matrix(kdf[,2:(ncol(kdf)-1)]),cluster_rows=FALSE,cluster_columns=TRUE,row_dend_reorder = TRUE,show_row_dend = FALSE)
h2<-Heatmap(as.matrix(kdf[,'hk_clusters']),cluster_columns=TRUE,cluster_rows=FALSE)

hlist <- h1+h2
tiff('PBandRA_RNAseq/output_figures/Hmap_PBandRA_clusters.tiff',width=800,height=2400,res=300)
draw(hlist, row_split = kdf$hk_clusters, cluster_row_slices = FALSE,cluster_rows=TRUE, row_dend_reorder = TRUE,show_row_dend = FALSE)    
dev.off()

#format and write output table
kdf_table<-kdf
kdf_table$gene<-rownames(kdf)
kdf_table<-kdf_table[,c('gene','hk_clusters')]
names(kdf_table)<-c('gene','cluster')
write.table(kdf_table,'PBandRA_RNAseq/output_data/PBandRA_clusters.txt',row.names=FALSE,quote=FALSE,col.names=TRUE)
       
       
# conduct gene ontology analysis of clusters ------------------------------------------

bg_universe<-allDE_collated$gene
length(bg_universe)

k<-5
clusters<-read.table('PBandRA_RNAseq/output_data/PBandRA_cluster_k5_complete.txt',stringsAsFactors = FALSE,header=TRUE)

clusters$cluster<-factor(clusters$cluster,levels=c('5','4','3','2','1'))
ck_go_cc<-compareCluster(gene~cluster, data=clusters[,c('gene','cluster')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='CC',pAdjustMethod='fdr', pvalueCutoff=0.01,universe=bg_universe,keyType='SYMBOL',minGSSize=20,maxGSSize=500)

dcc<-dotplot(ck_go_cc, showCategory = 10)+theme_bw()+theme(panel.grid.major.y = element_blank(),legend.position='none',plot.margin=unit(c(0,3,0,0),unit='cm'),
                                                           axis.text.x = element_text(angle = 315,vjust=0.5,hjust=0))+
  scale_colour_gradient(low='forestgreen',high='gold2')+coord_flip()+scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
ggsave(filename='PBandRA_RNAseq/output_figures/compareCluster_CC.png',dcc,width=12,height=4)

                                                                                      
# look at genes associated with differential patient survival -------------

 # data from R2 Genomics
NBsurvival_kscan<-read.table('PBandRA_RNAseq/data/survival_kscan.txt',stringsAsFactors = FALSE)
nrow(NBsurvival_kscan)
nrow(subset(NBsurvival_kscan,V5=='high'))
nrow(subset(NBsurvival_kscan,V5=='low'))
high<-subset(NBsurvival_kscan,V5=='high')$V2 #high expression associated with poor survival
low<-subset(NBsurvival_kscan,V5=='low')$V2 #low expression associated with poor survival

allDE_collated$survival_type<-rep('ns',nrow(allDE_collated))
for(i in 1:nrow(allDE_collated)){
  if((allDE_collated[i,'gene'] %in% high)==TRUE){
    allDE_collated[i,'survival_type']<-'high'
  }
  if((allDE_collated[i,'gene'] %in% low)==TRUE){
    allDE_collated[i,'survival_type']<-'low'
  }
}
allDE_collated$survival_type<-factor(allDE_collated$survival_type,levels=c('ns','low','high'))
allDE_collated<-allDE_collated[order(allDE_collated$survival_type),]

# volcano plot for PB+RA v DMSO RNA-seq data with survival data groups mapped on                                                                                      
png('PBandRA_RNAseq/output_figures/volcano_suvival_high.png',width=500,height=350)
ggplot(data=allDE_collated, aes(x=PBandRAvControl_log2FoldChange, y=-log10(PBandRAvControl_padj), col=survival_type))+geom_point()+theme_classic()+
  scale_color_manual(values=c("grey68","grey68","darkorange"))+geom_vline(xintercept=-0.5, col="black",linetype='dashed')+
  geom_hline(yintercept=-log10(0.05), col="black",linetype='dashed')+geom_vline(xintercept=0.5, col="black",linetype='dashed')+
  coord_trans(xlim=c(-6,6))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=20))+xlab(expression('log'[2]*'FC PB+RA v DMSO'))+ylab(expression("-log"[10]*"FDR"))
dev.off()

allDE_collated$survival_type<-factor(allDE_collated$survival_type,levels=c('ns','high','low'))
allDE_collated<-allDE_collated[order(allDE_collated$survival_type),]

png('PBandRA_RNAseq/output_figures/volcano_suvival_low.png',width=500,height=350)
ggplot(data=allDE_collated, aes(x=PBandRAvControl_log2FoldChange, y=-log10(PBandRAvControl_padj), col=survival_type))+geom_point()+theme_classic()+
  scale_color_manual(values=c("grey68","grey68","green4"))+geom_vline(xintercept=-0.5, col="black",linetype='dashed')+
  geom_hline(yintercept=-log10(0.05), col="black",linetype='dashed')+geom_vline(xintercept=0.5, col="black",linetype='dashed')+
  coord_trans(xlim=c(-6,6))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=20))+xlab(expression('log'[2]*'FC PB+RA v DMSO'))+ylab(expression("-log"[10]*"FDR"))
dev.off()


# statistics
dunnTest(PBandRAvControl_log2FoldChange~survival_type,data=allDE_collated,two.sided=TRUE)


