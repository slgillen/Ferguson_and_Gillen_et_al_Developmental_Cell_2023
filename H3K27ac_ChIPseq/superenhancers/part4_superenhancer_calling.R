#this script contains elements of code adapted from ROSE call_superenhancers.R


# load libraries ----------------------------------------------------------
library(ggplot2)



# import files ------------------------------------------------------------

#assumes run from within the directory containing consensus enhancer region maps

filenames <- list.files(pattern = "*consensus_ENHANCER_REGION_MAP.txt")
length(filenames)

BE2C_files<-filenames[grepl('BE2C',filenames)]
length(BE2C_files)
BE2C_files

IMR32_files<-filenames[grepl('IMR32',filenames)]
length(IMR32_files)
IMR32_files

SY5Y_files<-filenames[grepl('SY5Y',filenames)]
length(SY5Y_files)
SY5Y_files



# get average H3K27ac signal across the five biological replicates --------

############# BE2C ############
#collate data to get averages
BE2Cdata<-NULL
for(i in 1:length(BE2C_files)){
  Bdata<-read.delim(BE2C_files[i],stringsAsFactors = FALSE,header=TRUE)
  if(i<=1){
    BE2Cdata<-Bdata[,-5]
  }else{
    BE2Cdata<-merge(BE2Cdata,Bdata[,c(1,7)],by='REGION_ID')
  }
  rm(Bdata)
}

#get replicate averages
BE2Cdata$average_control<-apply(BE2Cdata[,grepl('control',names(BE2Cdata))==TRUE],1,function(x) mean(x,na.rm=TRUE))
BE2Cdata$average_PB<-apply(BE2Cdata[,grepl('PB',names(BE2Cdata))==TRUE],1,function(x) mean(x,na.rm=TRUE))


############# IMR32 ############
#collate data to get averages
IMR32data<-NULL
for(i in 1:length(IMR32_files)){
  Bdata<-read.delim(IMR32_files[i],stringsAsFactors = FALSE,header=TRUE)
  if(i<=1){
    IMR32data<-Bdata[,-5]
  }else{
    IMR32data<-merge(IMR32data,Bdata[,c(1,7)],by='REGION_ID')
  }
  rm(Bdata)
}

#get replicate averages
IMR32data$average_control<-apply(IMR32data[,grepl('control',names(IMR32data))==TRUE],1,function(x) mean(x,na.rm=TRUE))
IMR32data$average_PB<-apply(IMR32data[,grepl('PB',names(IMR32data))==TRUE],1,function(x) mean(x,na.rm=TRUE))


############# SH-SY5Y ############
#collate data to get averages
SY5Ydata<-NULL
for(i in 1:length(SY5Y_files)){
  Bdata<-read.delim(SY5Y_files[i],stringsAsFactors = FALSE,header=TRUE)
  if(i<=1){
    SY5Ydata<-Bdata[,-5]
  }else{
    SY5Ydata<-merge(SY5Ydata,Bdata[,c(1,7)],by='REGION_ID')
  }
  rm(Bdata)
}

#get replicate averages
SY5Ydata$average_control<-apply(SY5Ydata[,grepl('control',names(SY5Ydata))==TRUE],1,function(x) mean(x,na.rm=TRUE))
SY5Ydata$average_PB<-apply(SY5Ydata[,grepl('PB',names(SY5Ydata))==TRUE],1,function(x) mean(x,na.rm=TRUE))



# write output tables -----------------------------------------------------
BE2Cdata<-BE2Cdata[order(BE2Cdata$rank_control),]
write.table(BE2Cdata,'BE2C_supertable_full.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

IMR32data<-IMR32data[order(IMR32data$rank_control),]
write.table(IMR32data,'IMR32_supertable_full.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

SY5Ydata<-SY5Ydata[order(SY5Ydata$rank_control),]
write.table(SY5Ydata,'SY5Y_supertable_full.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')




# functions from ROSE script ---------------------------------------------------------------

#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#============================================================================

 #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
 calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}






# run ROSE script from each condition in each cell line --------------------------------------------------------------

############# BE2C ############
#average_control CPM
rankBy_vector <- as.numeric(BE2Cdata$average_control)

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab='',ylab='',lwd=2,col=4)
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)

plotFileName = paste('BE2C_control_RPM_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste('BE2C_control_RPM','_enhancers'),ylab=paste('H3K27ac',' Signal'),pch=19,cex=2)+
abline(h=cutoff_options$absolute,col='grey',lty=2)+abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)+lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')+
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)
dev.off()

rm(rankBy_vector, cutoff_options,superEnhancerRows)


#average_PB CPM
rankBy_vector <- as.numeric(BE2Cdata$average_PB)

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab='',ylab='',lwd=2,col=4)
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)

plotFileName = paste('BE2C_PB_RPM_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste('BE2C_PB_RPM','_enhancers'),ylab=paste('H3K27ac',' Signal'),pch=19,cex=2)+
abline(h=cutoff_options$absolute,col='grey',lty=2)+abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)+lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')+
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)
dev.off()

rm(rankBy_vector, cutoff_options,superEnhancerRows)



############# IMR32 ############
#average_control CPM
rankBy_vector <- as.numeric(IMR32data$average_control)

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab='',ylab='',lwd=2,col=4)
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)

plotFileName = paste('IMR32_control_RPM_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste('IMR32_control_RPM','_enhancers'),ylab=paste('H3K27ac',' Signal'),pch=19,cex=2)+
abline(h=cutoff_options$absolute,col='grey',lty=2)+abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)+lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')+
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)
dev.off()

rm(rankBy_vector, cutoff_options,superEnhancerRows)

#average_PB CPM
rankBy_vector <- as.numeric(IMR32data$average_PB)

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab='',ylab='',lwd=2,col=4)
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)

plotFileName = paste('IMR32_PB_RPM_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste('IMR32_PB_RPM','_enhancers'),ylab=paste('H3K27ac',' Signal'),pch=19,cex=2)+
abline(h=cutoff_options$absolute,col='grey',lty=2)+abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)+lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')+
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)
dev.off()

rm(rankBy_vector, cutoff_options,superEnhancerRows)


############# SY5Y ############
#average_control CPM
rankBy_vector <- as.numeric(SY5Ydata$average_control)

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab='',ylab='',lwd=2,col=4)
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)

plotFileName = paste('SY5Y_control_RPM_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste('SY5Y_control_RPM','_enhancers'),ylab=paste('H3K27ac',' Signal'),pch=19,cex=2)+
abline(h=cutoff_options$absolute,col='grey',lty=2)+abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)+lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')+
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)
dev.off()

rm(rankBy_vector, cutoff_options,superEnhancerRows)

#average_PB CPM
rankBy_vector <- as.numeric(SY5Ydata$average_PB)

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab='',ylab='',lwd=2,col=4)
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)

plotFileName = paste('SY5Y_PB_RPM_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste('SY5Y_PB_RPM','_enhancers'),ylab=paste('H3K27ac',' Signal'),pch=19,cex=2)+
abline(h=cutoff_options$absolute,col='grey',lty=2)+abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)+lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')+
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)
dev.off()

rm(rankBy_vector, cutoff_options,superEnhancerRows)




# re-call superenhancers across all cell lines and conditions --------------------------------------------------------------------
#ROSE initial script used to call super-enhancers, but for comparison between cell lines and conditions want a consistent threshold of signal to be used
#taken the minimum threshold found in the above analysis (60) and used this as a consistent threshold for superenhancer calling across conditions and cell lines

max(BE2Cdata$average_control)
max(BE2Cdata$average_PB)

max(IMR32data$average_control)
max(IMR32data$average_PB)

max(SY5Ydata$average_control)
max(SY5Ydata$average_PB)


#### BE2C ####
print(nrow(BE2Cdata))
BE2Cdata<-BE2Cdata[order(BE2Cdata$average_control,decreasing=TRUE),]
BE2Cdata$enhancerRank<-seq(1,nrow(BE2Cdata),1)
BE2Cdata$isSuper<-numeric(nrow(BE2Cdata))
BE2Cdata[BE2Cdata$average_control>60,'isSuper']<-1
BE2Cdata_ss<-subset(BE2Cdata,BE2Cdata[,'average_control']>60)
print(nrow(BE2Cdata_ss))

f1<-ggplot(BE2Cdata,aes(x=enhancerRank,y=average_control,col=factor(isSuper)))+geom_point()+theme_classic()+
scale_colour_manual(values=c('grey68','red3'))+ylab('H3K27ac signal (CPM)')+xlab('enhancer rank')+xlim(nrow(BE2Cdata),0)+
theme(legend.position = "none",axis.text.y=element_text(size=12),axis.text.x=element_text(size=10),axis.title=element_text(size=12))+coord_trans(ylim=c(0,2100))
ggsave(plot=f1,filename='BE2C_control_CPM_hockeystick.png',width=3.25,height=2.75)

BE2Cdata<-BE2Cdata[order(BE2Cdata$average_PB,decreasing=TRUE),]
BE2Cdata$enhancerRank<-seq(1,nrow(BE2Cdata),1)
BE2Cdata$isSuper<-numeric(nrow(BE2Cdata))
BE2Cdata[BE2Cdata$average_PB>60,'isSuper']<-1
BE2Cdata_ss<-subset(BE2Cdata,BE2Cdata[,'average_PB']>60)
print(nrow(BE2Cdata_ss))

f1<-ggplot(BE2Cdata,aes(x=enhancerRank,y=average_PB,col=factor(isSuper)))+geom_point()+theme_classic()+
scale_colour_manual(values=c('grey68','red3'))+ylab('H3K27ac signal (CPM)')+xlab('enhancer rank')+xlim(nrow(BE2Cdata),0)+
theme(legend.position = "none",axis.text.y=element_text(size=12),axis.text.x=element_text(size=10),axis.title=element_text(size=12))+coord_trans(ylim=c(0,2100))
ggsave(plot=f1,filename='BE2C_PB_CPM_hockeystick.png',width=3.25,height=2.75)


#### IMR32 ####
print(nrow(IMR32data))
IMR32data<-IMR32data[order(IMR32data$average_control,decreasing=TRUE),]
IMR32data$enhancerRank<-seq(1,nrow(IMR32data),1)
IMR32data$isSuper<-numeric(nrow(IMR32data))
IMR32data[IMR32data$average_control>60,'isSuper']<-1
IMR32data_ss<-subset(IMR32data,IMR32data[,'average_control']>60)
print(nrow(IMR32data_ss))

f1<-ggplot(IMR32data,aes(x=enhancerRank,y=average_control,col=factor(isSuper)))+geom_point()+theme_classic()+
scale_colour_manual(values=c('grey68','red3'))+ylab('H3K27ac signal (CPM)')+xlab('enhancer rank')+xlim(nrow(IMR32data),0)+
theme(legend.position = "none",axis.text.y=element_text(size=12),axis.text.x=element_text(size=10),axis.title=element_text(size=12))+coord_trans(ylim=c(0,1250))
ggsave(plot=f1,filename='IMR32_control_CPM_hockeystick.png',width=3.25,height=2.75)

IMR32data<-IMR32data[order(IMR32data$average_PB,decreasing=TRUE),]
IMR32data$enhancerRank<-seq(1,nrow(IMR32data),1)
IMR32data$isSuper<-numeric(nrow(IMR32data))
IMR32data[IMR32data$average_PB>60,'isSuper']<-1
IMR32data_ss<-subset(IMR32data,IMR32data[,'average_PB']>60)
print(nrow(IMR32data_ss))

f1<-ggplot(IMR32data,aes(x=enhancerRank,y=average_PB,col=factor(isSuper)))+geom_point()+theme_classic()+
scale_colour_manual(values=c('grey68','red3'))+ylab('H3K27ac signal (CPM)')+xlab('enhancer rank')+xlim(nrow(IMR32data),0)+
theme(legend.position = "none",axis.text.y=element_text(size=12),axis.text.x=element_text(size=10),axis.title=element_text(size=12))+coord_trans(ylim=c(0,1250))
ggsave(plot=f1,filename='IMR32_PB_CPM_hockeystick.png',width=3.25,height=2.75)


#### SY5Y ####
print(nrow(SY5Ydata))
SY5Ydata<-SY5Ydata[order(SY5Ydata$average_control,decreasing=TRUE),]
SY5Ydata$enhancerRank<-seq(1,nrow(SY5Ydata),1)
SY5Ydata$isSuper<-numeric(nrow(SY5Ydata))
SY5Ydata[SY5Ydata$average_control>60,'isSuper']<-1
SY5Ydata_ss<-subset(SY5Ydata,SY5Ydata[,'average_control']>60)
print(nrow(SY5Ydata_ss))

f1<-ggplot(SY5Ydata,aes(x=enhancerRank,y=average_control,col=factor(isSuper)))+geom_point()+theme_classic()+
scale_colour_manual(values=c('grey68','red3'))+ylab('H3K27ac signal (CPM)')+xlab('enhancer rank')+xlim(nrow(SY5Ydata),0)+
theme(legend.position = "none",axis.text.y=element_text(size=12),axis.text.x=element_text(size=10),axis.title=element_text(size=12))+coord_trans(ylim=c(0,3400))
ggsave(plot=f1,filename='SY5Y_control_CPM_hockeystick.png',width=3.25,height=2.75)

SY5Ydata<-SY5Ydata[order(SY5Ydata$average_PB,decreasing=TRUE),]
SY5Ydata$enhancerRank<-seq(1,nrow(SY5Ydata),1)
SY5Ydata$isSuper<-numeric(nrow(SY5Ydata))
SY5Ydata[SY5Ydata$average_PB>60,'isSuper']<-1
SY5Ydata_ss<-subset(SY5Ydata,SY5Ydata[,'average_PB']>60)
print(nrow(SY5Ydata_ss))

f1<-ggplot(SY5Ydata,aes(x=enhancerRank,y=average_PB,col=factor(isSuper)))+geom_point()+theme_classic()+
scale_colour_manual(values=c('grey68','red3'))+ylab('H3K27ac signal (CPM)')+xlab('enhancer rank')+xlim(nrow(SY5Ydata),0)+
theme(legend.position = "none",axis.text.y=element_text(size=12),axis.text.x=element_text(size=10),axis.title=element_text(size=12))+coord_trans(ylim=c(0,3400))
ggsave(plot=f1,filename='SY5Y_PB_CPM_hockeystick.png',width=3.25,height=2.75)







# classify superenhancers -------------------------------------------------

#classify superenhancers as increased, sustained or decreased


#### BE2C ####

BE2C_supertable<-read.table('BE2C_supertable_full.txt',stringsAsFactors = FALSE,header=TRUE)
names(BE2C_supertable)
dim(BE2C_supertable)

BE2C_supertable$control_SE<-0
BE2C_supertable[BE2C_supertable$average_control>60,'control_SE']<-1
BE2C_supertable$PB_SE<-0
BE2C_supertable[BE2C_supertable$average_PB>60,'PB_SE']<-1

BE2C_supertable$SE_type<-'none'
for(i in 1:nrow(BE2C_supertable)){
  if(BE2C_supertable[i,'control_SE']==1){
    if(BE2C_supertable[i,'PB_SE']==1){
      BE2C_supertable[i,'SE_type']<-'both'
    }else{
      BE2C_supertable[i,'SE_type']<-'control'
    }
  }else{
    if(BE2C_supertable[i,'PB_SE']==1){
      BE2C_supertable[i,'SE_type']<-'PB'
    }
  }
}
names(BE2Cdata)
BE2Cdata<-subset(BE2Cdata,SE_type!='none')
nrow(BE2Cdata) 

#get differential signal per replicate
BE2Cdata$Rep1_PBtoCon<-log(BE2Cdata$BE2C_Rep1_7dPB/BE2Cdata$BE2C_Rep1_control,2)
BE2Cdata$Rep2_PBtoCon<-log(BE2Cdata$BE2C_Rep2_7dPB/BE2Cdata$BE2C_Rep2_control,2)
BE2Cdata$Rep3_PBtoCon<-log(BE2Cdata$BE2C_Rep3_7dPB/BE2Cdata$BE2C_Rep3_control,2)
BE2Cdata$Rep4_PBtoCon<-log(BE2Cdata$BE2C_Rep4_7dPB/BE2Cdata$BE2C_Rep4_control,2)
BE2Cdata$Rep5_PBtoCon<-log(BE2Cdata$BE2C_Rep5_7dPB/BE2Cdata$BE2C_Rep5_control,2)

#get average signal differential
BE2Cdata$average_differential<-apply(BE2Cdata[,grepl('PBtoCon',names(BE2Cdata))==TRUE],1,function(x) mean(x,na.rm=TRUE))
summary(BE2Cdata$average_differential)
nrow(BE2Cdata)

#classify superenhancers
BE2Cdata$diff_type<-rep('none',nrow(BE2Cdata))
for(i in 1:nrow(BE2Cdata)){
  if(BE2Cdata[i,'average_differential']>0.15){
    BE2Cdata[i,'diff_type']<-'up'
  }else{
    if(BE2Cdata[i,'average_differential']<(-0.15)){
      BE2Cdata[i,'diff_type']<-'down'
    }else{
      BE2Cdata[i,'diff_type']<-'ns'
    }
  }
}

#create bed files for input into plotProfile
write.table(subset(BE2Cdata[,c(2,3,4,1,38,39)],diff_type=='up'),'BE2C_SE_up.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(subset(BE2Cdata[,c(2,3,4,1,38,39)],diff_type=='down'),'BE2C_SE_down.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(subset(BE2Cdata[,c(2,3,4,1,38,39)],diff_type=='ns'),'BE2C_SE_ns.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)


#create table for input into ROSE gene mapper
BE2C_supertableX<-BE2Cdata[,c(2,3,4,1,38,37,5,16,17,32:36)]
names(BE2C_supertableX)
supertable_names<-c('CHROM','START','STOP','REGION_ID','strand','NUM_LOCI','CONSTITUENT_SIZE','average_control','average_PB','average_differential','enhancerRank_control','enhancerRank_PB','enhancerRank_differential','IS_SUPER')
names(BE2C_supertableX)<-supertable_names
names(BE2C_supertableX)

BE2C_supertableX<-BE2C_supertableX[order(BE2C_supertableX$enhancerRank_differential),]
write.table(BE2C_supertableX,'BE2C_supertable.bed',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)




#### IMR32 ####

IMR32_supertable<-read.table('IMR32_supertable_full.txt',stringsAsFactors = FALSE,header=TRUE)
names(IMR32_supertable)
dim(IMR32_supertable)

IMR32_supertable$control_SE<-0
IMR32_supertable[IMR32_supertable$average_control>60,'control_SE']<-1
IMR32_supertable$PB_SE<-0
IMR32_supertable[IMR32_supertable$average_PB>60,'PB_SE']<-1

IMR32_supertable$SE_type<-'none'
for(i in 1:nrow(IMR32_supertable)){
  if(IMR32_supertable[i,'control_SE']==1){
    if(IMR32_supertable[i,'PB_SE']==1){
      IMR32_supertable[i,'SE_type']<-'both'
    }else{
      IMR32_supertable[i,'SE_type']<-'control'
    }
  }else{
    if(IMR32_supertable[i,'PB_SE']==1){
      IMR32_supertable[i,'SE_type']<-'PB'
    }
  }
}
names(IMR32data)
IMR32data<-subset(IMR32data,SE_type!='none')
nrow(IMR32data) 

#get differential signal per replicate
IMR32data$Rep1_PBtoCon<-log(IMR32data$IMR32_Rep1_5dPB/IMR32data$IMR32_Rep1_control,2)
IMR32data$Rep2_PBtoCon<-log(IMR32data$IMR32_Rep2_5dPB/IMR32data$IMR32_Rep2_control,2)
IMR32data$Rep3_PBtoCon<-log(IMR32data$IMR32_Rep3_5dPB/IMR32data$IMR32_Rep3_control,2)
IMR32data$Rep4_PBtoCon<-log(IMR32data$IMR32_Rep4_5dPB/IMR32data$IMR32_Rep4_control,2)
IMR32data$Rep5_PBtoCon<-log(IMR32data$IMR32_Rep5_5dPB/IMR32data$IMR32_Rep5_control,2)

#get average signal differential
IMR32data$average_differential<-apply(IMR32data[,grepl('PBtoCon',names(IMR32data))==TRUE],1,function(x) mean(x,na.rm=TRUE))
summary(IMR32data$average_differential)
nrow(IMR32data)

#classify superenhancers
IMR32data$diff_type<-rep('none',nrow(IMR32data))
for(i in 1:nrow(IMR32data)){
  if(IMR32data[i,'average_differential']>0.15){
    IMR32data[i,'diff_type']<-'up'
  }else{
    if(IMR32data[i,'average_differential']<(-0.15)){
      IMR32data[i,'diff_type']<-'down'
    }else{
      IMR32data[i,'diff_type']<-'ns'
    }
  }
}

#create bed files for input into plotProfile
write.table(subset(IMR32data[,c(2,3,4,1,38,39)],diff_type=='up'),'IMR32_SE_up.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(subset(IMR32data[,c(2,3,4,1,38,39)],diff_type=='down'),'IMR32_SE_down.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(subset(IMR32data[,c(2,3,4,1,38,39)],diff_type=='ns'),'IMR32_SE_ns.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)


#create table for input into ROSE gene mapper
IMR32_supertableX<-IMR32data[,c(2,3,4,1,38,37,5,16,17,32:36)]
names(IMR32_supertableX)
supertable_names<-c('CHROM','START','STOP','REGION_ID','strand','NUM_LOCI','CONSTITUENT_SIZE','average_control','average_PB','average_differential','enhancerRank_control','enhancerRank_PB','enhancerRank_differential','IS_SUPER')
names(IMR32_supertableX)<-supertable_names
names(IMR32_supertableX)

IMR32_supertableX<-IMR32_supertableX[order(IMR32_supertableX$enhancerRank_differential),]
write.table(IMR32_supertableX,'IMR32_supertable.bed',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)





#### SY5Y ####

SY5Y_supertable<-read.table('SY5Y_supertable_full.txt',stringsAsFactors = FALSE,header=TRUE)
names(SY5Y_supertable)
dim(SY5Y_supertable)

SY5Y_supertable$control_SE<-0
SY5Y_supertable[SY5Y_supertable$average_control>60,'control_SE']<-1
SY5Y_supertable$PB_SE<-0
SY5Y_supertable[SY5Y_supertable$average_PB>60,'PB_SE']<-1

SY5Y_supertable$SE_type<-'none'
for(i in 1:nrow(SY5Y_supertable)){
  if(SY5Y_supertable[i,'control_SE']==1){
    if(SY5Y_supertable[i,'PB_SE']==1){
      SY5Y_supertable[i,'SE_type']<-'both'
    }else{
      SY5Y_supertable[i,'SE_type']<-'control'
    }
  }else{
    if(SY5Y_supertable[i,'PB_SE']==1){
      SY5Y_supertable[i,'SE_type']<-'PB'
    }
  }
}
names(SY5Ydata)
SY5Ydata<-subset(SY5Ydata,SE_type!='none')
nrow(SY5Ydata) 

#get differential signal per replicate
SY5Ydata$Rep1_PBtoCon<-log(SY5Ydata$SY5Y_Rep1_5dPB/SY5Ydata$SY5Y_Rep1_control,2)
SY5Ydata$Rep2_PBtoCon<-log(SY5Ydata$SY5Y_Rep2_5dPB/SY5Ydata$SY5Y_Rep2_control,2)
SY5Ydata$Rep3_PBtoCon<-log(SY5Ydata$SY5Y_Rep3_5dPB/SY5Ydata$SY5Y_Rep3_control,2)
SY5Ydata$Rep4_PBtoCon<-log(SY5Ydata$SY5Y_Rep4_5dPB/SY5Ydata$SY5Y_Rep4_control,2)
SY5Ydata$Rep5_PBtoCon<-log(SY5Ydata$SY5Y_Rep5_5dPB/SY5Ydata$SY5Y_Rep5_control,2)

#get average signal differential
SY5Ydata$average_differential<-apply(SY5Ydata[,grepl('PBtoCon',names(SY5Ydata))==TRUE],1,function(x) mean(x,na.rm=TRUE))
summary(SY5Ydata$average_differential)
nrow(SY5Ydata)

#classify superenhancers
SY5Ydata$diff_type<-rep('none',nrow(SY5Ydata))
for(i in 1:nrow(SY5Ydata)){
  if(SY5Ydata[i,'average_differential']>0.15){
    SY5Ydata[i,'diff_type']<-'up'
  }else{
    if(SY5Ydata[i,'average_differential']<(-0.15)){
      SY5Ydata[i,'diff_type']<-'down'
    }else{
      SY5Ydata[i,'diff_type']<-'ns'
    }
  }
}

#create bed files for input into plotProfile
write.table(subset(SY5Ydata[,c(2,3,4,1,38,39)],diff_type=='up'),'SY5Y_SE_up.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(subset(SY5Ydata[,c(2,3,4,1,38,39)],diff_type=='down'),'SY5Y_SE_down.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
write.table(subset(SY5Ydata[,c(2,3,4,1,38,39)],diff_type=='ns'),'SY5Y_SE_ns.bed',quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)


#create table for input into ROSE gene mapper
SY5Y_supertableX<-SY5Ydata[,c(2,3,4,1,38,37,5,16,17,32:36)]
names(SY5Y_supertableX)
supertable_names<-c('CHROM','START','STOP','REGION_ID','strand','NUM_LOCI','CONSTITUENT_SIZE','average_control','average_PB','average_differential','enhancerRank_control','enhancerRank_PB','enhancerRank_differential','IS_SUPER')
names(SY5Y_supertableX)<-supertable_names
names(SY5Y_supertableX)

SY5Y_supertableX<-SY5Y_supertableX[order(SY5Y_supertableX$enhancerRank_differential),]
write.table(SY5Y_supertableX,'SY5Y_supertable.bed',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)

