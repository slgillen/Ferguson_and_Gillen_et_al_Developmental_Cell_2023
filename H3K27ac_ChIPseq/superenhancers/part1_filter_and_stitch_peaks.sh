#!/bin/bash




############ filter out the known amplified regions that are present in one or more of these cell lines ############ 
#MYCN, ALK and MEIS1

replicates='Rep1 Rep2 Rep3 Rep4 Rep5'

#SH-SY5Y
for rep in $replicates
do
  bedtools intersect -v -a $peakdir/SHSY5Y_control_${rep}_peaks.broadPeak -b $ampdir/amplified_regions_forexclusion.bed > $datadir/SY5Y_${rep}_control_peaks_filt.broadPeak &
done
wait

for rep in $replicates
do
  bedtools intersect -v -a $peakdir/SHSY5Y_5dPB_${rep}_peaks.broadPeak -b $ampdir/amplified_regions_forexclusion.bed > $datadir/SY5Y_${rep}_5dPB_peaks_filt.broadPeak &
done
wait

#IMR-32
for rep in $replicates
do
  bedtools intersect -v -a $indir/IMR32_control_${rep}_peaks.broadPeak -b $ampdir/amplified_regions_forexclusion.bed > $outdir/IMR32_${rep}_control_peaks_filt.broadPeak &
done
wait

for rep in $replicates
do
  bedtools intersect -v -a $indir/IMR32_5dPB_${rep}_peaks.broadPeak -b $ampdir/amplified_regions_forexclusion.bed > $outdir/IMR32_${rep}_5dPB_peaks_filt.broadPeak &
done
wait

#SK-N-BE(2)C
for rep in $replicates
do
  bedtools intersect -v -a $indir/BE2C_control_${rep}_peaks.broadPeak -b $ampdir/amplified_regions_forexclusion.bed > $outdir/BE2C_${rep}_control_peaks_filt.broadPeak &
done
wait

for rep in $replicates
do
  bedtools intersect -v -a $indir/BE2C_7dPB_${rep}_peaks.broadPeak -b $ampdir/amplified_regions_forexclusion.bed > $outdir/BE2C_${rep}_7dPB_peaks_filt.broadPeak &
done
wait


############################### format the data for region stitching ############################### 

#SH-SY5Y
for rep in $replicates
do
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $5, $2, $3, $7, $6, $9, $4}' $datadir/SY5Y_${rep}_control_peaks_filt.broadPeak > $datadir/SY5Y_${rep}_control_peaks_broad.gff &
done
wait

for rep in $replicates
do
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $5, $2, $3, $7, $6, $9, $4}' $datadir/SY5Y_${rep}_5dPB_peaks_filt.broadPeak > $datadir/SY5Y_${rep}_5dPB_peaks_broad.gff &
done
wait

#IMR-32
for rep in $replicates
do
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $5, $2, $3, $7, $6, $9, $4}' $outdir/IMR32_${rep}_control_peaks_filt.broadPeak > $outdir/IMR32_${rep}_control_peaks_broad.gff &
done
wait

for rep in $replicates
do
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $5, $2, $3, $7, $6, $9, $4}' $outdir/IMR32_${rep}_5dPB_peaks_filt.broadPeak > $outdir/IMR32_${rep}_5dPB_peaks_broad.gff &
done
wait


#SK-N-BE(2)C
for rep in $replicates
do
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $5, $2, $3, $7, $6, $9, $4}' $outdir/BE2C_${rep}_control_peaks_filt.broadPeak > $outdir/BE2C_${rep}_control_peaks_broad.gff &
done
wait

for rep in $replicates
do
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $4, $5, $2, $3, $7, $6, $9, $4}' $outdir/BE2C_${rep}_7dPB_peaks_filt.broadPeak > $outdir/BE2C_${rep}_7dPB_peaks_broad.gff &
done
wait

################################## stitch the broad peaks in to clusters using ROSE ################################## 
#max stitch distance 12.5kb; exclude TSS region 2.5kb

#SH-SY5Y
for rep in $replicates
do
  python stitching.py -g HG19 -i $datadir/SY5Y_${rep}_control_peaks_broad.gff -o $maindir/ -t 2500 -s 12500 &
done
wait

for rep in $replicates
do
  python stitching.py -g HG19 -i $datadir/SY5Y_${rep}_5dPB_peaks_broad.gff -o $maindir/ -t 2500 -s 12500 & 
done
wait


#IMR-32
for rep in $replicates
do
  python stitching.py -g HG19 -i $indir/IMR32_${rep}_control_peaks_broad.gff -o $maindir/ -t 2500 -s 12500 &
done
wait

for rep in $replicates
do
  python stitching.py -g HG19 -i $indir/IMR32_${rep}_5dPB_peaks_broad.gff -o $maindir/ -t 2500 -s 12500 & 
done
wait


#SK-N-BE(2)C
for rep in $replicates
do
  python stitching.py -g HG19 -i $indir/BE2C_${rep}_control_peaks_broad.gff -o $maindir/ -t 2500 -s 12500 & 
done
wait

for rep in $replicates
do
  python stitching.py -g HG19 -i $indir/BE2C_${rep}_7dPB_peaks_broad.gff -o $maindir/ -t 2500 -s 12500 &
done
wait
