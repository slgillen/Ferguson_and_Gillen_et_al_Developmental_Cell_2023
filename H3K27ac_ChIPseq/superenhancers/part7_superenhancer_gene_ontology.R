#this script conducts gene ontology on the genes that have been mapped in terms of proximity to the superenhancer regions using ROSE


# load libraries ----------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)
library(ggplot2)


# create output directory -------------------------------------------------
dir.create("PB_H3K27ac_ChIPseq/output_figures/SE_gene_onotology")


##########################################################################
############################# SK-N-BE(2)C ################################
##########################################################################

# import data  ------------------------------------------

#### BE2C ####
BE2Ctable<-read.table('BE2C_supertable.bed',stringsAsFactors = FALSE,header=FALSE)
supertable_names<-c('CHROM','START','STOP','REGION_ID','strand','NUM_LOCI','CONSTITUENT_SIZE','average_control','average_PB','average_differential','enhancerRank_control','enhancerRank_PB','enhancerRank_differential','IS_SUPER')
names(BE2Ctable)<-supertable_names
names(BE2Ctable)
summary(BE2Ctable$average_differential)
BE2Ctable$diff_type<-rep('none',nrow(BE2Ctable))
for(i in 1:nrow(BE2Ctable)){
  if(BE2Ctable[i,'average_differential']>0.15){
    BE2Ctable[i,'diff_type']<-'up'
  }else{
    if(BE2Ctable[i,'average_differential']<(-0.15)){
      BE2Ctable[i,'diff_type']<-'down'
    }else{
      if(BE2Ctable[i,'average_differential']<0.15 & BE2Ctable[i,'average_differential']>(-0.15)){
        BE2Ctable[i,'diff_type']<-'ns'
      }  
    }
  }
}



# import gene mappings ----------------------------------------------------
BE2C_GtoE<-read.delim('BE2C_supertable_GENE_TO_ENHANCER.txt',stringsAsFactors = FALSE,header=TRUE)

nrow(BE2C_GtoE) 
BE2C_GtoE_new<-NULL
for(i in 1:nrow(BE2C_GtoE)){
  nline<-BE2C_GtoE[i,]
  nenhancers<-strsplit(BE2C_GtoE[i,'PROXIMAL_ENHANCERS'],split=',')[[1]]
  if(length(nenhancers)==1){
    BE2C_GtoE_new<-rbind(BE2C_GtoE_new,nline)
  }else{
    for(j in 1:length(nenhancers)){
      nline2<-nline
      nline2[3]<-nenhancers[j]
      BE2C_GtoE_new<-rbind(BE2C_GtoE_new,nline2)
      nline2
    }
  }
}
names(BE2C_GtoE_new)<-names(BE2C_GtoE)
names(BE2C_GtoE_new)
BE2C_GtoE<-BE2C_GtoE_new

nrow(BE2C_GtoE)
length(unique(BE2C_GtoE$GENE_NAME))
BE2Ctable$ID1<-paste(BE2Ctable$CHROM,BE2Ctable$START,sep=':')
BE2Ctable$ID2<-paste(BE2Ctable$ID1,BE2Ctable$STOP,sep='-')
names(BE2Ctable)[17]<-'PROXIMAL_ENHANCERS'


# combine the tables ------------------------------------------------------
BE2C_GtoE<-full_join(BE2C_GtoE,BE2Ctable[,c('PROXIMAL_ENHANCERS','diff_type')],by='PROXIMAL_ENHANCERS')
nrow(BE2C_GtoE)
BE2C_GtoE<-subset(BE2C_GtoE,is.na(GENE_NAME)==FALSE)
nrow(BE2C_GtoE)


# get the gene sets ------------------------------------------------------

BE2C_up<-unique(subset(BE2C_GtoE,diff_type=='up')$GENE_NAME)
BE2C_ns<-unique(subset(BE2C_GtoE,diff_type=='ns')$GENE_NAME)
BE2C_down<-unique(subset(BE2C_GtoE,diff_type=='down')$GENE_NAME)

#make table for input
setnames<-c('up','ns','down')
BE2Cset<-list(BE2C_up,BE2C_ns,BE2C_down)
BE2Cdf<-NULL
for(i in 1:length(BE2Cset)){
  BE2Cdfx<-data.frame(BE2Cset[[i]],stringsAsFactors = FALSE)
  BE2Cdfx$diff_type<-rep(setnames[i],nrow(BE2Cdfx))
  BE2Cdf<-rbind(BE2Cdf,BE2Cdfx)
  rm(BE2Cdfx)
}
names(BE2Cdf)<-c('gene','diff_type')
BE2Cdf$diff_type<-factor(BE2Cdf$diff_type,levels=c('down','ns','up'))




# run the gene ontology analysis ---------------------------------------------------
bg_universe<-read.table('PB_H3K27ac_ChIPseq/data/GO_background.txt',stringsAsFactors = FALSE,header=FALSE) 
bg_universe<-bg_universe$V1


#plot BE2C gene ontology biological process and cellular component
go_all<-compareCluster(gene~diff_type, data=BE2Cdf[,c('gene','diff_type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='ALL',pAdjustMethod='fdr',universe=bg_universe, pvalueCutoff=0.1,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
dim(go_all@compareClusterResult)
go_all@compareClusterResult<-subset(go_all@compareClusterResult,ONTOLOGY!='MF')
dim(go_all@compareClusterResult)
go_all@compareClusterResult$diff_type<-factor(go_all@compareClusterResult$diff_type,levels=c('up','ns','down'))
go_all@compareClusterResult<-go_all@compareClusterResult[order(go_all@compareClusterResult$diff_type),]

dcc<-dotplot(go_all, x='diff_type',showCategory = 10)+theme_bw()+theme(panel.grid.major.x = element_blank())+scale_size(limits = c(0, 0.1))+
  scale_colour_gradient(low='forestgreen',high='gold2',limits=c(0,0.1))+scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  theme(text=element_text(family='sans'))
ggsave(filename='PB_H3K27ac_ChIPseq/output_figures/SE_gene_onotology/BE2C_SEs_GO.png',plot=dcc,width=5.5,height=4)
rm(go_all,dcc)



#########################################################################
############################### IMR-32 ##################################
#########################################################################

# import data  ------------------------------------------

#### IMR32 ####
IMR32table<-read.table('IMR32_supertable.bed',stringsAsFactors = FALSE,header=FALSE)
supertable_names<-c('CHROM','START','STOP','REGION_ID','strand','NUM_LOCI','CONSTITUENT_SIZE','average_control','average_PB','average_differential','enhancerRank_control','enhancerRank_PB','enhancerRank_differential','IS_SUPER')
names(IMR32table)<-supertable_names
names(IMR32table)
summary(IMR32table$average_differential)
IMR32table$diff_type<-rep('none',nrow(IMR32table))
for(i in 1:nrow(IMR32table)){
  if(IMR32table[i,'average_differential']>0.15){
    IMR32table[i,'diff_type']<-'up'
  }else{
    if(IMR32table[i,'average_differential']<(-0.15)){
      IMR32table[i,'diff_type']<-'down'
    }else{
      if(IMR32table[i,'average_differential']<0.15 & IMR32table[i,'average_differential']>(-0.15)){
        IMR32table[i,'diff_type']<-'ns'
      }  
    }
  }
}



# import gene mappings ----------------------------------------------------
IMR32_GtoE<-read.delim('IMR32_supertable_GENE_TO_ENHANCER.txt',stringsAsFactors = FALSE,header=TRUE)

nrow(IMR32_GtoE) 
IMR32_GtoE_new<-NULL
for(i in 1:nrow(IMR32_GtoE)){
  nline<-IMR32_GtoE[i,]
  nenhancers<-strsplit(IMR32_GtoE[i,'PROXIMAL_ENHANCERS'],split=',')[[1]]
  if(length(nenhancers)==1){
    IMR32_GtoE_new<-rbind(IMR32_GtoE_new,nline)
  }else{
    for(j in 1:length(nenhancers)){
      nline2<-nline
      nline2[3]<-nenhancers[j]
      IMR32_GtoE_new<-rbind(IMR32_GtoE_new,nline2)
      nline2
    }
  }
}
names(IMR32_GtoE_new)<-names(IMR32_GtoE)
names(IMR32_GtoE_new)
IMR32_GtoE<-IMR32_GtoE_new

nrow(IMR32_GtoE)
length(unique(IMR32_GtoE$GENE_NAME))
IMR32table$ID1<-paste(IMR32table$CHROM,IMR32table$START,sep=':')
IMR32table$ID2<-paste(IMR32table$ID1,IMR32table$STOP,sep='-')
names(IMR32table)[17]<-'PROXIMAL_ENHANCERS'


# combine the tables ------------------------------------------------------
IMR32_GtoE<-full_join(IMR32_GtoE,IMR32table[,c('PROXIMAL_ENHANCERS','diff_type')],by='PROXIMAL_ENHANCERS')
nrow(IMR32_GtoE)
IMR32_GtoE<-subset(IMR32_GtoE,is.na(GENE_NAME)==FALSE)
nrow(IMR32_GtoE)


# get the gene sets ------------------------------------------------------

IMR32_up<-unique(subset(IMR32_GtoE,diff_type=='up')$GENE_NAME)
IMR32_ns<-unique(subset(IMR32_GtoE,diff_type=='ns')$GENE_NAME)
IMR32_down<-unique(subset(IMR32_GtoE,diff_type=='down')$GENE_NAME)

#make table for input
setnames<-c('up','ns','down')
IMR32set<-list(IMR32_up,IMR32_ns,IMR32_down)
IMR32df<-NULL
for(i in 1:length(IMR32set)){
  IMR32dfx<-data.frame(IMR32set[[i]],stringsAsFactors = FALSE)
  IMR32dfx$diff_type<-rep(setnames[i],nrow(IMR32dfx))
  IMR32df<-rbind(IMR32df,IMR32dfx)
  rm(IMR32dfx)
}
names(IMR32df)<-c('gene','diff_type')
IMR32df$diff_type<-factor(IMR32df$diff_type,levels=c('down','ns','up'))




# run the gene ontology analysis ---------------------------------------------------
bg_universe<-read.table('PB_H3K27ac_ChIPseq/data/GO_background.txt',stringsAsFactors = FALSE,header=FALSE) 
bg_universe<-bg_universe$V1


#plot IMR32 gene ontology biological process and cellular component
go_all<-compareCluster(gene~diff_type, data=IMR32df[,c('gene','diff_type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='ALL',pAdjustMethod='fdr',universe=bg_universe, pvalueCutoff=0.1,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
dim(go_all@compareClusterResult)
go_all@compareClusterResult<-subset(go_all@compareClusterResult,ONTOLOGY!='MF')
dim(go_all@compareClusterResult)
go_all@compareClusterResult$diff_type<-factor(go_all@compareClusterResult$diff_type,levels=c('up','ns','down'))
go_all@compareClusterResult<-go_all@compareClusterResult[order(go_all@compareClusterResult$diff_type),]

dcc<-dotplot(go_all, x='diff_type',showCategory = 10)+theme_bw()+theme(panel.grid.major.x = element_blank())+scale_size(limits = c(0, 0.1))+
  scale_colour_gradient(low='forestgreen',high='gold2',limits=c(0,0.1))+scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  theme(text=element_text(family='sans'))
ggsave(filename='PB_H3K27ac_ChIPseq/output_figures/SE_gene_onotology/IMR32_SEs_GO.png',plot=dcc,width=5.5,height=4)
rm(go_all,dcc)




#########################################################################
############################### SH-SY5Y #################################
#########################################################################

# import data  ------------------------------------------

#### SY5Y ####
SY5Ytable<-read.table('SY5Y_supertable.bed',stringsAsFactors = FALSE,header=FALSE)
supertable_names<-c('CHROM','START','STOP','REGION_ID','strand','NUM_LOCI','CONSTITUENT_SIZE','average_control','average_PB','average_differential','enhancerRank_control','enhancerRank_PB','enhancerRank_differential','IS_SUPER')
names(SY5Ytable)<-supertable_names
names(SY5Ytable)
summary(SY5Ytable$average_differential)
SY5Ytable$diff_type<-rep('none',nrow(SY5Ytable))
for(i in 1:nrow(SY5Ytable)){
  if(SY5Ytable[i,'average_differential']>0.15){
    SY5Ytable[i,'diff_type']<-'up'
  }else{
    if(SY5Ytable[i,'average_differential']<(-0.15)){
      SY5Ytable[i,'diff_type']<-'down'
    }else{
      if(SY5Ytable[i,'average_differential']<0.15 & SY5Ytable[i,'average_differential']>(-0.15)){
        SY5Ytable[i,'diff_type']<-'ns'
      }  
    }
  }
}



# import gene mappings ----------------------------------------------------
SY5Y_GtoE<-read.delim('SY5Y_supertable_GENE_TO_ENHANCER.txt',stringsAsFactors = FALSE,header=TRUE)

nrow(SY5Y_GtoE) 
SY5Y_GtoE_new<-NULL
for(i in 1:nrow(SY5Y_GtoE)){
  nline<-SY5Y_GtoE[i,]
  nenhancers<-strsplit(SY5Y_GtoE[i,'PROXIMAL_ENHANCERS'],split=',')[[1]]
  if(length(nenhancers)==1){
    SY5Y_GtoE_new<-rbind(SY5Y_GtoE_new,nline)
  }else{
    for(j in 1:length(nenhancers)){
      nline2<-nline
      nline2[3]<-nenhancers[j]
      SY5Y_GtoE_new<-rbind(SY5Y_GtoE_new,nline2)
      nline2
    }
  }
}
names(SY5Y_GtoE_new)<-names(SY5Y_GtoE)
names(SY5Y_GtoE_new)
SY5Y_GtoE<-SY5Y_GtoE_new

nrow(SY5Y_GtoE)
length(unique(SY5Y_GtoE$GENE_NAME))
SY5Ytable$ID1<-paste(SY5Ytable$CHROM,SY5Ytable$START,sep=':')
SY5Ytable$ID2<-paste(SY5Ytable$ID1,SY5Ytable$STOP,sep='-')
names(SY5Ytable)[17]<-'PROXIMAL_ENHANCERS'


# combine the tables ------------------------------------------------------
SY5Y_GtoE<-full_join(SY5Y_GtoE,SY5Ytable[,c('PROXIMAL_ENHANCERS','diff_type')],by='PROXIMAL_ENHANCERS')
nrow(SY5Y_GtoE)
SY5Y_GtoE<-subset(SY5Y_GtoE,is.na(GENE_NAME)==FALSE)
nrow(SY5Y_GtoE)


# get the gene sets ------------------------------------------------------

SY5Y_up<-unique(subset(SY5Y_GtoE,diff_type=='up')$GENE_NAME)
SY5Y_ns<-unique(subset(SY5Y_GtoE,diff_type=='ns')$GENE_NAME)
SY5Y_down<-unique(subset(SY5Y_GtoE,diff_type=='down')$GENE_NAME)

#make table for input
setnames<-c('up','ns','down')
SY5Yset<-list(SY5Y_up,SY5Y_ns,SY5Y_down)
SY5Ydf<-NULL
for(i in 1:length(SY5Yset)){
  SY5Ydfx<-data.frame(SY5Yset[[i]],stringsAsFactors = FALSE)
  SY5Ydfx$diff_type<-rep(setnames[i],nrow(SY5Ydfx))
  SY5Ydf<-rbind(SY5Ydf,SY5Ydfx)
  rm(SY5Ydfx)
}
names(SY5Ydf)<-c('gene','diff_type')
SY5Ydf$diff_type<-factor(SY5Ydf$diff_type,levels=c('down','ns','up'))




# run the gene ontology analysis ---------------------------------------------------
bg_universe<-read.table('PB_H3K27ac_ChIPseq/data/GO_background.txt',stringsAsFactors = FALSE,header=FALSE) 
bg_universe<-bg_universe$V1


#plot SY5Y gene ontology biological process and cellular component
go_all<-compareCluster(gene~diff_type, data=SY5Ydf[,c('gene','diff_type')], fun="enrichGO", OrgDb="org.Hs.eg.db", ont='ALL',pAdjustMethod='fdr',universe=bg_universe, pvalueCutoff=0.1,keyType='SYMBOL',minGSSize=20,maxGSSize=500)
dim(go_all@compareClusterResult)
go_all@compareClusterResult<-subset(go_all@compareClusterResult,ONTOLOGY!='MF')
dim(go_all@compareClusterResult)
go_all@compareClusterResult$diff_type<-factor(go_all@compareClusterResult$diff_type,levels=c('up','ns','down'))
go_all@compareClusterResult<-go_all@compareClusterResult[order(go_all@compareClusterResult$diff_type),]

dcc<-dotplot(go_all, x='diff_type',showCategory = 10)+theme_bw()+theme(panel.grid.major.x = element_blank())+scale_size(limits = c(0, 0.1))+
  scale_colour_gradient(low='forestgreen',high='gold2',limits=c(0,0.1))+scale_y_discrete(labels = function(x) str_wrap(x, width = 100))+
  theme(text=element_text(family='sans'))
ggsave(filename='PB_H3K27ac_ChIPseq/output_figures/SE_gene_onotology/SY5Y_SEs_GO.png',plot=dcc,width=5.5,height=4)
rm(go_all,dcc)

