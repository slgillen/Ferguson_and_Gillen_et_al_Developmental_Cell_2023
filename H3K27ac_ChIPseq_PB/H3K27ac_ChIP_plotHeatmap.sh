#!/bin/bash

#plot regions of H3K27ac marks determined to be differential by DiffBind analysis
#non-sigificantly changed regions also included for comparison

#SK-N-BE(2)C

computeMatrix reference-point --referencePoint center -b 2500 -a 2500 -R $datadir/BE2C_broadpeak_up05.bed $datadir/BE2C_broadpeak_ns05.bed $datadir/BE2C_broadpeak_down05.bed -S $bwdir/BE2C_control_H3K27ac_ChIP_mean_CPM.bw $bwdir/BE2C_7dPB_H3K27ac_ChIP_mean_CPM.bw -p 10 \
-o $datadir/computeMatrix_BE2C05.gz --outFileSortedRegions $datadir/computeMatrix_BE2C05.bed --sortRegions no


plotHeatmap -m $datadir/computeMatrix_BE2C05.gz -out $datadir/hmap_BE2C05_x2.png \
 --colorMap PuRd PuRd --whatToShow 'heatmap and colorbar' --sortRegions no --missingDataColor=white \
 --xAxisLabel '' --heatmapHeight 17 --heatmapWidth 3 --zMin 0 0 --zMax 1.6 1.6



#IMR-32

computeMatrix reference-point --referencePoint center -b 2500 -a 2500 -R $datadir/IMR32_broadpeak_up05.bed $datadir/IMR32_broadpeak_ns05.bed $datadir/IMR32_broadpeak_down05.bed -S $bwdir/IMR32_control_H3K27ac_ChIP_mean_CPM.bw $bwdir/IMR32_5dPB_H3K27ac_ChIP_mean_CPM.bw -p 10 \
-o $datadir/computeMatrix_IMR3205.gz --outFileSortedRegions $datadir/computeMatrix_IMR3205.bed --sortRegions no

plotHeatmap -m $datadir/computeMatrix_IMR3205.gz -out $datadir/hmap_IMR3205_x2.png \
 --colorMap YlOrBr YlOrBr --whatToShow 'heatmap and colorbar' --sortRegions no --missingDataColor=white \
 --xAxisLabel '' --heatmapHeight 17 --heatmapWidth 3 --legendLocation lower-center --zMin 0 0 --zMax 1.51 1.51



#SH-SY5Y

computeMatrix reference-point --referencePoint center -b 2500 -a 2500 -R $datadir/SY5Y_broadpeak_up05.bed $datadir/SY5Y_broadpeak_ns05.bed $datadir/SY5Y_broadpeak_down05.bed -S $bwdir/SY5Y_control_H3K27ac_ChIP_mean_CPM.bw $bwdir/SY5Y_5dPB_H3K27ac_ChIP_mean_CPM.bw -p 10 \
-o $datadir/computeMatrix_SY5Y05.gz --outFileSortedRegions $datadir/computeMatrix_SY5Y05.bed --sortRegions no

plotHeatmap -m $datadir/computeMatrix_SY5Y05.gz -out $datadir/hmap_SY5Y05_x2.png \
 --colorMap PuBu PuBu --whatToShow 'heatmap and colorbar' --sortRegions no --missingDataColor=white \
 --xAxisLabel '' --heatmapHeight 17 --heatmapWidth 3 --legendLocation lower-center --zMin 0 0 --zMax 1.9 1.9

