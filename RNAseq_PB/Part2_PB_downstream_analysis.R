#this script conducts downstream investigations of the palbociclib RNA-seq data in three neuroblastoma cell lines


# load libraries ----------------------------------------------------------
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(factoextra)
library(scales)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(fgsea)


# collate datasets --------------------------------------------------------

# combine differential expression tables from the different comparisons conducted
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


data_path='PB_RNAseq/output_data/'
filenames <- dir(data_path, pattern = "*.csv")
length(filenames)

allinputs<-lapply(as.list(filenames),function(x) read.delim(paste0(data_path,x),sep=',',stringsAsFactors = FALSE,header=TRUE))
length(allinputs)

allDE_collated<-collate_results(allinputs, filenames)
names(allDE_collated)
dim(allDE_collated)

                  

# look at E2F target genes ------------------------------------------------

pathwaystouse <- gmtPathways('PB_RNAseq/data/h.all.v7.1.symbols.gmt')
E2Ftargets<-pathwaystouse$HALLMARK_E2F_TARGETS
length(E2Ftargets)

E2Fdata<-allDE_collated[,grepl('gene|24hPBvControl',names(allDE_collated))]
dim(E2Fdata)
names(E2Fdata)

E2Fdata$E2Ftargets<-rep('all other',nrow(E2Fdata))
for(i in 1:nrow(E2Fdata)){
  if((E2Fdata[i,'gene'] %in% E2Ftargets)==TRUE){
    E2Fdata[i,'E2Ftargets']<-'E2F target'
  }
}
nrow(subset(E2Fdata,E2Ftargets=='E2F target'))

E2Fdata$E2Ftargets<-factor(E2Fdata$E2Ftargets,levels=c('all other','E2F target'))
E2Fdata<-E2Fdata[order(E2Fdata$E2Ftargets),]

f1<-ggplot(E2Fdata,aes(x=E2Ftargets,fill=E2Ftargets,y=SHSY5Y_24hPBvControl_log2FoldChange))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('grey90','grey40'))+xlab('')+ylab(expression('log'[2]*'FC 24h PB v Control'))+
  theme(axis.text=element_text(size=14,colour='black'),axis.title=element_text(size=14,colour='black'),legend.position = 'none')+coord_trans(ylim=c(-5.5,3))
ggsave(plot=f1,filename='PB_RNAseq/output_figures/SY5Y_E2Ftargets.tiff',width=2.75,height=4,dpi=300)

f2<-ggplot(E2Fdata,aes(x=E2Ftargets,fill=E2Ftargets,y=BE2C_24hPBvControl_log2FoldChange))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('grey90','grey40'))+xlab('')+ylab(expression('log'[2]*'FC 24h PB v Control'))+
  theme(axis.text=element_text(size=14,colour='black'),axis.title=element_text(size=14,colour='black'),legend.position = 'none')+coord_trans(ylim=c(-5.5,3))
ggsave(plot=f2,filename='PB_RNAseq/output_figures/BE2C_E2Ftargets.tiff',width=2.75,height=4,dpi=300)

f3<-ggplot(E2Fdata,aes(x=E2Ftargets,fill=E2Ftargets,y=IMR32_24hPBvControl_log2FoldChange))+geom_boxplot(outlier.size=0.5)+theme_classic()+
  scale_fill_manual(values=c('grey90','grey40'))+xlab('')+ylab(expression('log'[2]*'FC 24h PB v Control'))+
  theme(axis.text=element_text(size=14,colour='black'),axis.title=element_text(size=14,colour='black'),legend.position = 'none')+coord_trans(ylim=c(-5.5,3))
ggsave(plot=f3,filename='PB_RNAseq/output_figures/IMR32_E2Ftargets.tiff',width=2.75,height=4,dpi=300)


             
# select differential genes ------------------------------------------------

unique_names<-names(allDE_collated)
unique_names<-unique(gsub('_log2FoldChange|_padj','',unique_names))
unique_names<-unique_names[unique_names!='gene']

for(n in 1:length(unique_names)){
  allDE_collated[,paste0(unique_names[n], '_SIGchange', collapse='')]<-rep('ns',nrow(allDE_collated))
}

#singificant change classified as p.adj<0.05 & log2FC > 0.5 or < -0.5
for(n in 1:length(unique_names)){
  for(i in 1:nrow(allDE_collated)){
    if(is.na(allDE_collated[i,paste0(unique_names[n], '_padj', collapse='')])==FALSE & allDE_collated[i,paste0(unique_names[n], '_padj', collapse='')]<0.05){
      if(allDE_collated[i,paste0(unique_names[n], '_log2FoldChange', collapse='')]>0.5){
        allDE_collated[i,paste0(unique_names[n], '_SIGchange', collapse='')]<-'up'
      }
      if(allDE_collated[i,paste0(unique_names[n], '_log2FoldChange', collapse='')]<(-0.5)){
        allDE_collated[i,paste0(unique_names[n], '_SIGchange', collapse='')]<-'down'
      }
    }
  }
}


IMR32up<-subset(allDE_collated,IMR32_24hPBvControl_SIGchange=='up'|IMR32_5dPBvControl_SIGchange=='up'|IMR32_5dPBv24hPB_SIGchange=='up')$gene
length(IMR32up)
IMR32down<-subset(allDE_collated,IMR32_24hPBvControl_SIGchange=='down'|IMR32_5dPBvControl_SIGchange=='down'|IMR32_5dPBv24hPB_SIGchange=='down')$gene
length(IMR32down)

BE2Cup<-subset(allDE_collated,BE2C_24hPBvControl_SIGchange=='up'|BE2C_7dPBvControl_SIGchange=='up'|BE2C_7dPBv24hPB_SIGchange=='up')$gene
length(BE2Cup)
BE2Cdown<-subset(allDE_collated,BE2C_24hPBvControl_SIGchange=='down'|BE2C_7dPBvControl_SIGchange=='down'|BE2C_7dPBv24hPB_SIGchange=='down')$gene
length(BE2Cdown)

SHSY5Yup<-subset(allDE_collated,SHSY5Y_24hPBvControl_SIGchange=='up'|SHSY5Y_5dPBvControl_SIGchange=='up'|SHSY5Y_5dPBv24hPB_SIGchange=='up')$gene
length(SHSY5Yup)
SHSY5Ydown<-subset(allDE_collated,SHSY5Y_24hPBvControl_SIGchange=='down'|SHSY5Y_5dPBvControl_SIGchange=='down'|SHSY5Y_5dPBv24hPB_SIGchange=='down')$gene
length(SHSY5Ydown)

overlap_up<-IMR32up[(IMR32up %in% BE2Cup)==TRUE]
overlap_up<-overlap_up[(overlap_up %in% SHSY5Yup)==TRUE]
length(overlap_up)

overlap_down<-IMR32down[(IMR32down %in% BE2Cdown)==TRUE]
overlap_down<-overlap_down[(overlap_down %in% SHSY5Ydown)==TRUE]
length(overlap_down)

all_overlap<-unique(c(overlap_up,overlap_down))
length(all_overlap)

         
# get CPM normalised values -----------------------------------------------

BE2C_RNAseq_counts<-read.csv('PB_RNAseq/data/BE2C_PB_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
IMR32_RNAseq_counts<-read.csv('PB_RNAseq/data/IMR32_PB_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
SHSY5Y_RNAseq_counts<-read.csv('PB_RNAseq/data/SHSY5Y_PB_RNAseq_rawcounts.csv',stringsAsFactors=FALSE,header=TRUE)
                  
             
BE2C_RNAseq_CPM<-BE2C_RNAseq_counts
BE2C_RNAseq_CPM[,2:ncol(BE2C_RNAseq_CPM)]<-apply(BE2C_RNAseq_CPM[,2:ncol(BE2C_RNAseq_CPM)], 2, function(x) (x/sum(x,na.rm=TRUE))*1000000)

IMR32_RNAseq_CPM<-IMR32_RNAseq_counts
IMR32_RNAseq_CPM[,2:ncol(IMR32_RNAseq_CPM)]<-apply(IMR32_RNAseq_CPM[,2:ncol(IMR32_RNAseq_CPM)], 2, function(x) (x/sum(x,na.rm=TRUE))*1000000)

SHSY5Y_RNAseq_CPM<-SHSY5Y_RNAseq_counts
SHSY5Y_RNAseq_CPM[,2:ncol(SHSY5Y_RNAseq_CPM)]<-apply(SHSY5Y_RNAseq_CPM[,2:ncol(SHSY5Y_RNAseq_CPM)], 2, function(x) (x/sum(x,na.rm=TRUE))*1000000)
                

# kmeans clustering -------------------------------------------------------

# filter for only differential gene set
BE2C_RNAseq_CPM<-subset(BE2C_RNAseq_CPM,(gene %in% all_overlap)==TRUE)
IMR32_RNAseq_CPM<-subset(IMR32_RNAseq_CPM,(gene %in% all_overlap)==TRUE)
SHSY5Y_RNAseq_CPM<-subset(SHSY5Y_RNAseq_CPM,(gene %in% all_overlap)==TRUE)
                              
# conduct scaling within each cell line
scale_df<-function(df_CPM){
  df_CPM_prescale<-df_CPM
  for(i in 1:nrow(df_CPM_prescale)){
    df_CPM_prescale[i,2:ncol(df_CPM_prescale)]<-scale(as.numeric(df_CPM_prescale[i,2:ncol(df_CPM_prescale)]))
  }
  return(df_CPM_prescale)
}

SHSY5Y_scaled<-scale_df(SHSY5Y_RNAseq_CPM)
IMR32_scaled<-scale_df(IMR32_RNAseq_CPM)
BE2C_scaled<-scale_df(BE2C_RNAseq_CPM)
         
# combined scaled data from all cell lines
scaled_merge<-merge(BE2C_scaled,IMR32_scaled,by='gene')
scaled_merge<-merge(scaled_merge,SHSY5Y_scaled,by='gene')
dim(scaled_merge)

# format data
kdf<-scaled_merge[,2:ncol(scaled_merge)]
rownames(kdf)<-scaled_merge$gene

# determine optimal number of clusters
#NOTE: with the selection of significantly differentially expressed either up or down regulated genes here, two clusters will dominate first
png('PB_RNAseq/output_figures/alllines_kmeans_clusters_elbowmethod.png',width=400,height=250)
print(fviz_nbclust(kdf, kmeans, method='wss'))
dev.off()

png('PB_RNAseq/output_figures/alllines_kmeans_clusters_elbowmethod_zoom.png',width=400,height=250)
print(fviz_nbclust(kdf, kmeans, method='wss'))+coord_trans(ylim=c(25000,40000))
dev.off()

# visual check of the value of k choice
k<-6
res.k<-hkmeans(kdf,k,hc.method='complete')
png(paste0('PB_RNAseq/output_figures/alllines_kmeans_',k,'_clusters.png'),width=350,height=300)
print(fviz_cluster(res.k, ellipse=FALSE,labelsize=6,ggtheme = theme_classic(),geom=c('point'),pointsize=1,stand=TRUE,shape='circle')+scale_colour_manual(values = c("darkorange", "green4","mediumvioletred",'royalblue3','darkgoldenrod2','red3','grey78','grey28')))
dev.off()
                  
# conduct clustering
res_k<-hkmeans(kdf,k,hc.method='complete')
kdf$k_clusters<-res_k$cluster
kdf$k_clusters<-as.character(kdf$k_clusters)
kdf<-kdf[order(kdf$hk_clusters),]

# extract number of genes in each cluster
for(i in 1:k){
  print(nrow(subset(kdf,hk_clusters==i)))
}                                                         
                  
# format and write output table
kdf_table<-kdf
kdf_table$gene<-rownames(kdf)
kdf_table<-kdf_table[,c('gene','hk_clusters')]
names(kdf_table)<-c('gene','cluster')
write.table(kdf_table,file=paste0('PB_RNAseq/output_data/cluster_details_',k,'_alllines.txt'),row.names=FALSE,quote=FALSE,col.names=TRUE)

             
                  
# plot clusters as heatmap ------------------------------------------------

h1<-Heatmap(as.matrix(kdf[,1:(ncol(kdf)-1)]),cluster_rows=TRUE,cluster_columns=FALSE,use_raster=FALSE,row_dend_reorder = TRUE,show_row_dend = FALSE,show_row_names=FALSE,
            column_split=factor(c(rep('BE2C',15),rep('IMR32',15),rep('SHSY5Y',15)),levels=c('BE2C','IMR32','SHSY5Y')),show_column_names = FALSE)
h2<-Heatmap(as.matrix(kdf[,'hk_clusters']),cluster_columns=FALSE,cluster_rows=FALSE)

hlist <- h1+h2
tiff(paste0('PB_RNAseq/output_figures/alllines_Hmap_kmeans',k,'_clusters.tiff'),width=1000,height=2600,res=300)
draw(hlist, row_split = kdf$hk_clusters, cluster_row_slices = TRUE,cluster_rows=FALSE, row_dend_reorder = TRUE,show_row_dend = FALSE)
dev.off()
               
                                                     
# gene ontology analysis of clusters --------------------------------------

bg_universe<-allDE_collated$gene
length(bg_universe)

clusters<-read.table('PB_RNAseq/output_data/cluster_details_6_alllines.txt',stringsAsFactors = FALSE,header=TRUE)

GOtypes<-c('CC')
for(GOtype in GOtypes){
  for(i in 1:k){
    geneList<-subset(clusters,cluster==i)$gene
    ego <- enrichGO(geneList, OrgDb = "org.Hs.eg.db", ont=GOtype, pAdjustMethod='fdr', pvalueCutoff=0.01,universe=bg_universe,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
    ego_df<-data.frame(ego,stringsAsFactors = FALSE)
    ego_df<-ego_df[order(ego_df$p.adjust),]
    write.csv(ego_df[,c(2,3,4,6,8)],file = paste0('clusters_',i,'_',GOtype,'.csv'),quote=FALSE,row.names=FALSE)
    ego_df$neglog10padj<--log(ego_df$p.adjust,10)
    ego_df<-ego_df[order(ego_df$p.adjust,decreasing=TRUE),] #order table by p-value
    print(nrow(ego_df))
    
    ego_df$Description<-factor(ego_df$Description,level=ego_df$Description)
    
    ego_df$gene_ratio<-sapply(ego_df$GeneRatio, function(x) as.numeric(strsplit(x,split='/')[[1]][1])/as.numeric(strsplit(x,split='/')[[1]][2]))
    
    #if more than 10 significant terms only plot the 10 most significant
    if(nrow(ego_df)>=10){
      png(paste0('clusters',k,'_',i,'_',GOtype,'.png'),width=450,height=175) 
      print(ggplot(ego_df[((nrow(ego_df)-9):nrow(ego_df)),],aes(x=neglog10padj,y=Description,fill=gene_ratio))+geom_bar(stat='identity')+theme_classic()+
              scale_x_continuous(expand=c(0,0))+labs(fill = 'Gene Ratio')+xlab('-log10 adjusted p-val')+ylab('')+theme(axis.text=element_text(size=12),axis.title=element_text(size=12)))
      dev.off()
    }else{
      png(paste0('clusters',k,'_',i,'_',GOtype,'.png'),width=450,height=175) 
      print(ggplot(ego_df,aes(x=neglog10padj,y=Description,fill=gene_ratio))+geom_bar(stat='identity')+theme_classic()+
              scale_x_continuous(expand=c(0,0))+labs(fill = 'Gene Ratio')+xlab('-log10 adjusted p-val')+ylab('')+theme(axis.text=element_text(size=12),axis.title=element_text(size=12)))
      dev.off()
    }
    rm(ego,ego_df,geneList)
  }
  
}
            
                                                     
                                                     
                                                     
                                                     
                                                     
                                                     
                                                     
                                                     
                                                     
                                                     

