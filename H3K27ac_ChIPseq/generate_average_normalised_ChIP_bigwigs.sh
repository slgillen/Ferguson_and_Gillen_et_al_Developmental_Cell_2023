#!/bin/bash

#this script takes the sorted bam files and create CPM normalised bigwig files
#then it get a mean normalised coverage for each condition (mean of the 5 biological replicates)

###########################################################################
############################### SK-N-BE(2)C ###############################
###########################################################################

#get all filename starts
filenames=$(find . -name BE2C\*\.sorted\.bam | sed 's/.\///' | sed s/\.sorted\.bam//)
echo $filenames

###### convert bam to CPM normalised bigwig for each sample ######
covdir=$maindir/bam_coverage
mkdir $covdir

for file in $filenames
do
  echo $file
  bamCoverage -b $indir/sorted_bam/${file}.sorted.bam -o $covdir/${file}_CPM.bw --binSize 50 --normalizeUsing CPM \
    --effectiveGenomeSize 2864785220 --extendReads &
done
wait


###### merge replicates ######
mergedir=$maindir/bam_coverage/merged_bw
mkdir $mergedir

### control samples ###

#get mean across bigwigs
wiggletools mean $covdir/BE2C*control_ChIP_CPM.bw > $mergedir/BE2C_control_CPM_merged.bedgraph

#sort the file
awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/BE2C_control_CPM_merged.bedgraph > $mergedir/BE2C_control_CPM_merged2.bedgraph
sort -k1,1 -k2,2n $mergedir/BE2C_control_CPM_merged2.bedgraph > $mergedir/BE2C_control_CPM_merged.sorted.bedgraph 

#convert to bigwig
bedGraphToBigWig $mergedir/BE2C_control_CPM_merged.sorted.bedgraph $covdir/chr_order.txt $mergedir/BE2C_control_H3K27ac_ChIP_mean_CPM.bw  

#remove unneeded intermediate files 
rm $mergedir/BE2C_control_CPM_merged.bedgraph
rm $mergedir/BE2C_control_CPM_merged2.bedgraph
rm $mergedir/BE2C_control_CPM_merged.sorted.bedgraph

### PB-treated samples ###
#get mean across bigwigs
wiggletools mean $covdir/BE2C*7dPB_ChIP_CPM.bw > $mergedir/BE2C_7dPB_CPM_merged.bedgraph

#sort the file
awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/BE2C_7dPB_CPM_merged.bedgraph > $mergedir/BE2C_7dPB_CPM_merged2.bedgraph
sort -k1,1 -k2,2n $mergedir/BE2C_7dPB_CPM_merged2.bedgraph > $mergedir/BE2C_7dPB_CPM_merged.sorted.bedgraph 

#convert to bigwig
bedGraphToBigWig $mergedir/BE2C_7dPB_CPM_merged.sorted.bedgraph $covdir/chr_order.txt $mergedir/BE2C_7dPB_H3K27ac_ChIP_mean_CPM.bw  

#remove unneeded intermediate files 
rm $mergedir/BE2C_7dPB_CPM_merged.bedgraph
rm $mergedir/BE2C_7dPB_CPM_merged2.bedgraph
rm $mergedir/BE2C_7dPB_CPM_merged.sorted.bedgraph


###### bigwig compare ######

bigwigCompare -b1 $mergedir/BE2C_7dPB_H3K27ac_ChIP_mean_CPM.bw -b2 $mergedir/BE2C_control_H3K27ac_ChIP_mean_CPM.bw --operation subtract -p 10 --outFileFormat bigwig \
 --outFileName $mergedir/BE2C_7dPBtoControl_subtract.bw




##########################################################################
################################# IMR-32 #################################
##########################################################################

#get all filename starts
filenames=$(find . -name IMR32\*\.sorted\.bam | sed 's/.\///' | sed s/\.sorted\.bam//)
echo $filenames

###### convert bam to CPM normalised bigwig for each sample ######
covdir=$maindir/bam_coverage
mkdir $covdir

for file in $filenames
do
  echo $file
  bamCoverage -b $indir/sorted_bam/${file}.sorted.bam -o $covdir/${file}_CPM.bw --binSize 50 --normalizeUsing CPM \
    --effectiveGenomeSize 2864785220 --extendReads &
done
wait

###### merge replicates ######
mergedir=$maindir/bam_coverage/merged_bw
mkdir $mergedir

### control samples ###

#get mean across bigwigs
wiggletools mean $covdir/IMR32*control_ChIP_CPM.bw > $mergedir/IMR32_control_CPM_merged.bedgraph

#sort the file
awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/IMR32_control_CPM_merged.bedgraph > $mergedir/IMR32_control_CPM_merged2.bedgraph
sort -k1,1 -k2,2n $mergedir/IMR32_control_CPM_merged2.bedgraph > $mergedir/IMR32_control_CPM_merged.sorted.bedgraph 

#convert to bigwig
bedGraphToBigWig $mergedir/IMR32_control_CPM_merged.sorted.bedgraph $covdir/chr_order.txt $mergedir/IMR32_control_H3K27ac_ChIP_mean_CPM.bw  

#remove unneeded intermediate files 
rm $mergedir/IMR32_control_CPM_merged.bedgraph
rm $mergedir/IMR32_control_CPM_merged2.bedgraph
rm $mergedir/IMR32_control_CPM_merged.sorted.bedgraph

### PB-treated samples ###
#get mean across bigwigs
wiggletools mean $covdir/IMR32*5dPB_ChIP_CPM.bw > $mergedir/IMR32_5dPB_CPM_merged.bedgraph

#sort the file
awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/IMR32_5dPB_CPM_merged.bedgraph > $mergedir/IMR32_5dPB_CPM_merged2.bedgraph
sort -k1,1 -k2,2n $mergedir/IMR32_5dPB_CPM_merged2.bedgraph > $mergedir/IMR32_5dPB_CPM_merged.sorted.bedgraph 

#convert to bigwig
bedGraphToBigWig $mergedir/IMR32_5dPB_CPM_merged.sorted.bedgraph $covdir/chr_order.txt $mergedir/IMR32_5dPB_H3K27ac_ChIP_mean_CPM.bw  

#remove unneeded intermediate files 
rm $mergedir/IMR32_5dPB_CPM_merged.bedgraph
rm $mergedir/IMR32_5dPB_CPM_merged2.bedgraph
rm $mergedir/IMR32_5dPB_CPM_merged.sorted.bedgraph


###### bigwig compare ######

bigwigCompare -b1 $mergedir/IMR32_5dPB_H3K27ac_ChIP_mean_CPM.bw -b2 $mergedir/IMR32_control_H3K27ac_ChIP_mean_CPM.bw --operation subtract -p 10 --outFileFormat bigwig \
 --outFileName $mergedir/IMR32_5dPBtoControl_subtract.bw


 

###########################################################################
################################# SH-SY5Y #################################
###########################################################################

#get all filename starts
filenames=$(find . -name SHSY5Y\*\.sorted\.bam | sed 's/.\///' | sed s/\.sorted\.bam//)
echo $filenames

###### convert bam to CPM normalised bigwig for each sample ######
covdir=$maindir/bam_coverage
mkdir $covdir

for file in $filenames
do
  echo $file
  bamCoverage -b $indir/sorted_bam/${file}.sorted.bam -o $covdir/${file}_CPM.bw --binSize 50 --normalizeUsing CPM \
    --effectiveGenomeSize 2864785220 --extendReads &
done
wait

###### merge replicates ######
mergedir=$maindir/bam_coverage/merged_bw
mkdir $mergedir

### control samples ###

#get mean across bigwigs
wiggletools mean $covdir/SHSY5Y*control_ChIP_CPM.bw > $mergedir/SHSY5Y_control_CPM_merged.bedgraph

#sort the file
awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/SHSY5Y_control_CPM_merged.bedgraph > $mergedir/SHSY5Y_control_CPM_merged2.bedgraph
sort -k1,1 -k2,2n $mergedir/SHSY5Y_control_CPM_merged2.bedgraph > $mergedir/SHSY5Y_control_CPM_merged.sorted.bedgraph 

#convert to bigwig
bedGraphToBigWig $mergedir/SHSY5Y_control_CPM_merged.sorted.bedgraph $covdir/chr_order.txt $mergedir/SHSY5Y_control_H3K27ac_ChIP_mean_CPM.bw  

#remove unneeded intermediate files 
rm $mergedir/SHSY5Y_control_CPM_merged.bedgraph
rm $mergedir/SHSY5Y_control_CPM_merged2.bedgraph
rm $mergedir/SHSY5Y_control_CPM_merged.sorted.bedgraph

### PB-treated samples ###
#get mean across bigwigs
wiggletools mean $covdir/SHSY5Y*5dPB_ChIP_CPM.bw > $mergedir/SHSY5Y_5dPB_CPM_merged.bedgraph

#sort the file
awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/SHSY5Y_5dPB_CPM_merged.bedgraph > $mergedir/SHSY5Y_5dPB_CPM_merged2.bedgraph
sort -k1,1 -k2,2n $mergedir/SHSY5Y_5dPB_CPM_merged2.bedgraph > $mergedir/SHSY5Y_5dPB_CPM_merged.sorted.bedgraph 

#convert to bigwig
bedGraphToBigWig $mergedir/SHSY5Y_5dPB_CPM_merged.sorted.bedgraph $covdir/chr_order.txt $mergedir/SHSY5Y_5dP_H3K27ac_ChIP_mean_CPM.bw  

#remove unneeded intermediate files 
rm $mergedir/SHSY5Y_5dPB_CPM_merged.bedgraph
rm $mergedir/SHSY5Y_5dPB_CPM_merged2.bedgraph
rm $mergedir/SHSY5Y_5dPB_CPM_merged.sorted.bedgraph


###### bigwig compare ######

bigwigCompare -b1 $mergedir/SHSY5Y_5dPB_H3K27ac_ChIP_mean_CPM.bw -b2 $mergedir/SHSY5Y_control_H3K27ac_ChIP_mean_CPM.bw --operation subtract -p 10 --outFileFormat bigwig \
 --outFileName $mergedir/SHSY5Y_5dPBtoControl_subtract.bw
