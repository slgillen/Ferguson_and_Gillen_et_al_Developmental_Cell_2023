# load libraries ----------------------------------------------------------

library(ggplot2)
library(clusterProfiler)

# combine the differential H3K27ac data from the different comparisons --------------------


# sort mark groups into non-redundant groups  ---------------------------

#special case antagonism groups
PBup_RAdown_PBRAup<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type=='up' & RAvDMSO_sig_type=='down')
PBup_RAdown_PBRAdown<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type=='down' & RAvDMSO_sig_type=='down')
PBup_RAdown_PBRAns<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type=='ns' & RAvDMSO_sig_type=='down')
PBdown_RAup_PBRAdown<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type=='down' & RAvDMSO_sig_type=='up')
PBdown_RAup_PBRAup<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type=='up' & RAvDMSO_sig_type=='up')
PBdown_RAup_PBRAns<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type=='ns' & RAvDMSO_sig_type=='up')

antag_regions<-c(PBup_RAdown_PBRAup$REGION_ID, PBup_RAdown_PBRAdown$REGION_ID, PBup_RAdown_PBRAns$REGION_ID, PBdown_RAup_PBRAdown$REGION_ID, PBdown_RAup_PBRAup$REGION_ID, PBdown_RAup_PBRAns$REGION_ID)

#up regulated regions
up_all3<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type=='up' & RAvDMSO_sig_type=='up' & (REGION_ID %in% antag_regions)==FALSE)
up_PB_PBRA<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type=='up' & RAvDMSO_sig_type!='up' & (REGION_ID %in% antag_regions)==FALSE)
up_RA_PBRA<-subset(mergedDB,PBvDMSO_sig_type!='up' & PBRAvDMSO_sig_type=='up' & RAvDMSO_sig_type=='up' & (REGION_ID %in% antag_regions)==FALSE)
up_PB_RA<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type!='up' & RAvDMSO_sig_type=='up' & (REGION_ID %in% antag_regions)==FALSE)
up_PB<-subset(mergedDB,PBvDMSO_sig_type=='up' & PBRAvDMSO_sig_type!='up' & RAvDMSO_sig_type!='up' & (REGION_ID %in% antag_regions)==FALSE)
up_PBRA<-subset(mergedDB,PBvDMSO_sig_type!='up' & PBRAvDMSO_sig_type=='up' & RAvDMSO_sig_type!='up' & (REGION_ID %in% antag_regions)==FALSE)
up_RA<-subset(mergedDB,PBvDMSO_sig_type!='up' & PBRAvDMSO_sig_type!='up' & RAvDMSO_sig_type=='up' & (REGION_ID %in% antag_regions)==FALSE)

#downregulated regions
down_all3<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type=='down' & RAvDMSO_sig_type=='down' & (REGION_ID %in% antag_regions)==FALSE)
down_PB_PBRA<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type=='down' & RAvDMSO_sig_type!='down' & (REGION_ID %in% antag_regions)==FALSE)
down_RA_PBRA<-subset(mergedDB,PBvDMSO_sig_type!='down' & PBRAvDMSO_sig_type=='down' & RAvDMSO_sig_type=='down' & (REGION_ID %in% antag_regions)==FALSE)
down_PB_RA<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type!='down' & RAvDMSO_sig_type=='down' & (REGION_ID %in% antag_regions)==FALSE)
down_PB<-subset(mergedDB,PBvDMSO_sig_type=='down' & PBRAvDMSO_sig_type!='down' & RAvDMSO_sig_type!='down' & (REGION_ID %in% antag_regions)==FALSE)
down_PBRA<-subset(mergedDB,PBvDMSO_sig_type!='down' & PBRAvDMSO_sig_type=='down' & RAvDMSO_sig_type!='down' & (REGION_ID %in% antag_regions)==FALSE)
down_RA<-subset(mergedDB,PBvDMSO_sig_type!='down' & PBRAvDMSO_sig_type!='down' & RAvDMSO_sig_type=='down' & (REGION_ID %in% antag_regions)==FALSE)

#regions with no change
ns_all3<-subset(mergedDB,PBvDMSO_sig_type=='ns' & PBRAvDMSO_sig_type=='ns' & RAvDMSO_sig_type=='ns' & (REGION_ID %in% antag_regions)==FALSE)

up_all3$DiffGroup<-rep('up_all3',nrow(up_all3))
up_PB_PBRA$DiffGroup<-rep('up_PB_PBRA',nrow(up_PB_PBRA))
up_RA_PBRA$DiffGroup<-rep('up_RA_PBRA',nrow(up_RA_PBRA))
up_PB_RA$DiffGroup<-rep('up_PB_RA',nrow(up_PB_RA))
up_PB$DiffGroup<-rep('up_PB',nrow(up_PB))
up_PBRA$DiffGroup<-rep('up_PBRA',nrow(up_PBRA))
up_RA$DiffGroup<-rep('up_RA',nrow(up_RA))

down_all3$DiffGroup<-rep('down_all3',nrow(down_all3))
down_PB_PBRA$DiffGroup<-rep('down_PB_PBRA',nrow(down_PB_PBRA))
down_RA_PBRA$DiffGroup<-rep('down_RA_PBRA',nrow(down_RA_PBRA))
down_PB_RA$DiffGroup<-rep('down_PB_RA',nrow(down_PB_RA))
down_PB$DiffGroup<-rep('down_PB',nrow(down_PB))
down_PBRA$DiffGroup<-rep('down_PBRA',nrow(down_PBRA))
down_RA$DiffGroup<-rep('down_RA',nrow(down_RA))

ns_all3$DiffGroup<-rep('ns_all3',nrow(ns_all3))

PBup_RAdown_PBRAup$DiffGroup<-rep('PBup_RAdown_PBRAup',nrow(PBup_RAdown_PBRAup))
PBup_RAdown_PBRAdown$DiffGroup<-rep('PBup_RAdown_PBRAdown',nrow(PBup_RAdown_PBRAdown))
PBup_RAdown_PBRAns$DiffGroup<-rep('PBup_RAdown_PBRAns',nrow(PBup_RAdown_PBRAns))
PBdown_RAup_PBRAdown$DiffGroup<-rep('PBdown_RAup_PBRAdown',nrow(PBdown_RAup_PBRAdown))
PBdown_RAup_PBRAup$DiffGroup<-rep('PBdown_RAup_PBRAup',nrow(PBdown_RAup_PBRAup))
PBdown_RAup_PBRAns$DiffGroup<-rep('PBdown_RAup_PBRAns',nrow(PBdown_RAup_PBRAns))

#combine all group types
allplotdata<-rbind(up_all3,up_PB_PBRA,up_RA_PBRA,up_PBRA, up_PB_RA,up_PB, up_RA, down_all3, down_PB_PBRA, down_RA_PBRA,down_PBRA, down_PB_RA, down_PB, down_RA, ns_all3,
                   PBup_RAdown_PBRAup,PBup_RAdown_PBRAdown,PBup_RAdown_PBRAns,PBdown_RAup_PBRAdown,PBdown_RAup_PBRAup,PBdown_RAup_PBRAns)

str(allplotdata)
allplotdata$DiffGroup<-factor(allplotdata$DiffGroup)
str(allplotdata)

# plot the H3K27ac mark group changes for each condition compared to DMSO ---------------------------------------------------------------

p1<-allplotdata[,c(1,3,21)]
names(p1)[2]<-'Fold'
p1$FoldType<-rep('PBvDMSO',nrow(p1))
p2<-allplotdata[,c(1,6,21)]
names(p2)[2]<-'Fold'
p2$FoldType<-rep('PBRAvDMSO',nrow(p2))
p3<-allplotdata[,c(1,9,21)]
names(p3)[2]<-'Fold'
p3$FoldType<-rep('RAvDMSO',nrow(p3))
p4<-allplotdata[,c(1,12,21)]
names(p4)[2]<-'Fold'
p4$FoldType<-rep('PBvRAvPB',nrow(p4))
p5<-allplotdata[,c(1,15,21)]
names(p5)[2]<-'Fold'
p5$FoldType<-rep('PBRAvRA',nrow(p5))
p6<-allplotdata[,c(1,18,21)]
names(p6)[2]<-'Fold'
p6$FoldType<-rep('PBvRA',nrow(p6))

allplotdata2<-rbind(p1,p2,p3,p4,p5,p6)
rm(p1,p2,p3,p4,p5,p6)
dim(allplotdata2)
allplotdata2$FoldType<-factor(allplotdata2$FoldType)
allplotdata2$DiffGroup<-factor(allplotdata2$DiffGroup,levels=c("down_all3","down_RA_PBRA","down_PB_PBRA","down_PBRA","down_RA","down_PB","down_PB_RA","ns_all3",
                                                               "PBdown_RAup_PBRAdown","PBdown_RAup_PBRAns","PBdown_RAup_PBRAup","PBup_RAdown_PBRAdown","PBup_RAdown_PBRAns","PBup_RAdown_PBRAup",
                                                               "up_PB_RA","up_PB","up_RA","up_PBRA","up_PB_PBRA","up_RA_PBRA","up_all3"))

outdir<-'final_region_plots/'

plotdf<-subset(allplotdata2,grepl('DMSO',allplotdata2$FoldType)==TRUE)
plotdf$FoldType<-factor(plotdf$FoldType,levels=c('PBvDMSO','RAvDMSO','PBRAvDMSO'))
f1<-ggplot(plotdf,aes(x=DiffGroup, y=Fold, fill=FoldType))+geom_boxplot(outlier.size=0.01,width=0.75)+theme_classic()+
  scale_fill_manual(values=c('red3','royalblue3','purple'))+theme(axis.text.x = element_text(angle = 315,vjust=0.5,hjust=0),legend.title=element_blank(),
                                                                  axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))+
  ylab('log2FC H3K27ac signal')+coord_trans(ylim=c(-3,4))
ggsave(filename=paste0(outdir,'H3K27ac_groups_boxplot.png'),plot=f1,width=13.5,height=4.5)
rm(plotdf,f1)


# connect H3K27ac regions to genes ---------------------------------------------------------------


# gene ontology on associated genes ---------------------------------------------------------------

GOdir<-'final_region_plots/gene_ontology/'

allplotdata<-merge(allplotdata,peakAnno[,c('REGION_ID','SYMBOL')],by='REGION_ID')
dim(allplotdata)

plotgroups<-c('up_PBRA','up_PB_PBRA','up_RA_PBRA','up_all3')
GOdata<-subset(allplotdata[,c('SYMBOL','DiffGroup')],(DiffGroup %in% plotgroups)==TRUE)
GOdata$DiffGroup<-factor(GOdata$DiffGroup,levels=plotgroups)

ck_go_cc<-compareCluster(SYMBOL~DiffGroup, data=GOdata, fun="enrichGO", OrgDb="org.Hs.eg.db", ont='CC',pAdjustMethod='fdr', pvalueCutoff=0.05,keyType='SYMBOL',minGSSize=10,maxGSSize=500)

dcc<-dotplot(ck_go_cc, showCategory = 10)+theme_bw()+theme(panel.grid.major.y = element_blank(),legend.position='none',plot.margin=unit(c(0,3,0,0),unit='cm'),
                                                           axis.text.x = element_text(angle = 315,vjust=0.5,hjust=0))+
  scale_colour_gradient(low='forestgreen',high='gold2')+coord_flip()+scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
ggsave(filename=paste0(GOdir,'/compareCluster_CC.png'),dcc,width=12,height=4)

rm(ck_go_cc)
rm(dcc)
















