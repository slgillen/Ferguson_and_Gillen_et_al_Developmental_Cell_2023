#!/bin/bash

#this script plots profiles of the average H3K27ac signal in the control and PB-treated condition for each cell line  

superdir=''
bwdir=''
profiledir=$maindir/SE_profiles
mkdir $profiledir

############################################## IMR-32 ##############################################
xname='IMR32_SEs'
computeMatrix scale-regions -R $superdir/IMR32_SE_up.bed $superdir/IMR32_SE_ns.bed $superdir/IMR32_SE_down.bed -S $bwdir/IMR32_control_H3K27ac_ChIP_mean_CPM.bw $bwdir/IMR32_5dPB_H3K27ac_ChIP_mean_CPM.bw \
-o $profiledir/computeMatrix_${xname}.gz --outFileSortedRegions $profiledir/computeMatrix_${xname}.bed -b 0 -a 0 --sortRegions keep


plotProfile -m $profiledir/computeMatrix_${xname}_nooutliers.gz -out $profiledir/profile_${xname}.png  \
  --plotTitle "IMR32 H3K27ac" --refPointLabel "center" -T "H3K27ac CPM" --colors grey \#DA8C13 red \
  --plotHeight 7 --plotWidth 7 --averageType mean --legendLocation lower-center --samplesLabel control PB --yMin 0 --yMax 1.4

plotProfile -m $profiledir/computeMatrix_${xname}_nooutliers.gz -out $profiledir/profile_${xname}_pergroup.png  \
  --perGroup --plotTitle "IMR32 H3K27ac" --refPointLabel "center" -T "H3K27ac CPM" --colors grey \#DA8C13 red \
  --plotHeight 6 --plotWidth 5.5 --averageType mean --legendLocation lower-center --samplesLabel control PB --yMin 0 --yMax 1.4


############################################## SK-N-BE(2)C ##############################################
xname='BE2C_SEs'
computeMatrix scale-regions -R $superdir/BE2C_SE_up.bed $superdir/BE2C_SE_ns.bed $superdir/BE2C_SE_down.bed -S $bwdir/BE2C_control_H3K27ac_ChIP_mean_CPM.bw $bwdir/BE2C_7dPB_H3K27ac_ChIP_mean_CPM.bw \
-o $profiledir/computeMatrix_${xname}.gz --outFileSortedRegions $profiledir/computeMatrix_${xname}.bed -b 0 -a 0 --sortRegions keep

plotProfile -m $profiledir/computeMatrix_${xname}_nooutliers.gz -out $profiledir/profile_${xname}.png  \
  --plotTitle "BE2C H3K27ac" --refPointLabel "center" -T "H3K27ac CPM" --colors grey \#B8155F red \
  --plotHeight 7 --plotWidth 7 --averageType mean --legendLocation lower-center --samplesLabel control PB --yMin 0 --yMax 1.1

plotProfile -m $profiledir/computeMatrix_${xname}_nooutliers.gz -out $profiledir/profile_${xname}_pergroup.png  \
  --perGroup --plotTitle "BE2C H3K27ac" --refPointLabel "center" -T "H3K27ac CPM" --colors grey \#B8155F red \
  --plotHeight 6 --plotWidth 5.5 --averageType mean --legendLocation lower-center --samplesLabel control PB --yMin 0 --yMax 1.1


############################################## SH-SY5Y ##############################################
xname='SY5Y_SEs'
computeMatrix scale-regions -R $superdir/SY5Y_SE_up.bed $superdir/SY5Y_SE_ns.bed $superdir/SY5Y_SE_down.bed -S $bwdir/SY5Y_control_H3K27ac_ChIP_mean_CPM.bw $bwdir/SY5Y_5dPB_H3K27ac_ChIP_mean_CPM.bw \
-o $profiledir/computeMatrix_${xname}.gz --outFileSortedRegions $profiledir/computeMatrix_${xname}.bed -b 0 -a 0 --sortRegions keep

plotProfile -m $profiledir/computeMatrix_${xname}_nooutliers.gz -out $profiledir/profile_${xname}.png  \
  --plotTitle "SY5Y H3K27ac" --refPointLabel "center" -T "H3K27ac CPM" --colors grey \#147FB1 red \
  --plotHeight 7 --plotWidth 7 --averageType mean --legendLocation lower-center --samplesLabel control PB --yMin 0 --yMax 1.5

plotProfile -m $profiledir/computeMatrix_${xname}_nooutliers.gz -out $profiledir/profile_${xname}_pergroup.png  \
  --perGroup --plotTitle "SY5Y H3K27ac" --refPointLabel "center" -T "H3K27ac CPM" --colors grey \#147FB1 red \
  --plotHeight 6 --plotWidth 5.5 --averageType mean --legendLocation lower-center --samplesLabel control PB --yMin 0 --yMax 1.5


