#!/bin/bash

#this script averages the bigwig files by a mean across the replicates for each condition

#using WiggleTools version 1.2


#################### directories #####################
maindir='' #path to main directory

filtdir=$maindir/bam_filtering
bwdir=$maindir/bigwig
mergedir=$maindir/merged_bigwig
mkdir $mergedir


###### conditions ######
#could set to be read as input from command line
#infuture ensure users use no extra _ between parts of sample name 
conditions='DMSO RA PB RA_PB'

############################### merge replicates (mean) ###########################
##### for IPs #####

for cond in ${conditions[@]}
do
    echo ${cond}
    wiggletools mean $bwdir/BE2C_${cond}_*_IP_CPM.bw > $mergedir/BE2C_mean_${cond}_IP_CPM.bedgraph
    awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/BE2C_mean_${cond}_IP_CPM.bedgraph > $mergedir/BE2C_mean_${cond}_IP_CPM_x2.bedgraph
    sort -k1,1 -k2,2n $mergedir/BE2C_mean_${cond}_IP_CPM_x2.bedgraph > $mergedir/BE2C_mean_${cond}_IP_CPM_sorted.bedgraph 
    bedGraphToBigWig $mergedir/BE2C_mean_${cond}_IP_CPM_sorted.bedgraph $maindir/core_files/UCSC_hg19_chromosome_sizes.txt $mergedir/BE2C_mean_${cond}_IP_CPM.bw 

    rm $mergedir/BE2C_mean_${cond}_IP_CPM.bedgraph
    rm $mergedir/BE2C_mean_${cond}_IP_CPM_x2.bedgraph
    rm $mergedir/BE2C_mean_${cond}_IP_CPM_sorted.bedgraph
done
wait

##### for inputs #####

for cond in ${conditions[@]}
do
    wiggletools mean $bwdir/BE2C_${cond}_*_Input_CPM.bw > $mergedir/BE2C_mean_${cond}_Input_CPM.bedgraph
    awk '(NF==4){print $0}' OFS='\t' FS='\t' $mergedir/BE2C_mean_${cond}_Input_CPM.bedgraph > $mergedir/BE2C_mean_${cond}_Input_CPM_x2.bedgraph
    sort -k1,1 -k2,2n $mergedir/BE2C_mean_${cond}_Input_CPM_x2.bedgraph > $mergedir/BE2C_mean_${cond}_Input_CPM_sorted.bedgraph 
    bedGraphToBigWig $mergedir/BE2C_mean_${cond}_Input_CPM_sorted.bedgraph $maindir/core_files/UCSC_hg19_chromosome_sizes.txt $mergedir/BE2C_mean_${cond}_Input_CPM.bw 

    rm $mergedir/BE2C_mean_${cond}_Input_CPM.bedgraph
    rm $mergedir/BE2C_mean_${cond}_Input_CPM_x2.bedgraph
    rm $mergedir/BE2C_mean_${cond}_Input_CPM_sorted.bedgraph
done
wait








