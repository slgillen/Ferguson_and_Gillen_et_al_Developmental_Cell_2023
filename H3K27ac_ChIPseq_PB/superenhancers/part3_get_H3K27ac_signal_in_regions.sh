#!/bin/bash

#this script formats the file attributes for use with featureCounts
#then uses featureCounts to get the counts of the paired end reads that overlap the consensus regions obtained in part 2
#the counts are then normalised to account for library size

maindir=''
indir=$maindir/merged_stitched_regions
outdir=$maindir/featureCounts
mkdir $outdir
logdir=$maindir/logs
mkdir $logdir


############################## sort attributes ##############################

python sort_attributes.py $indir/SY5Y_all_cat_btmerge.bed $indir/SY5Y_all_cat_btmerge.gff
python sort_attributes.py $indir/IMR32_all_cat_btmerge.bed $indir/IMR32_all_cat_btmerge.gff
python sort_attributes.py $indir/BE2C_all_cat_btmerge.bed $indir/BE2C_all_cat_btmerge.gff


######################### get paired-end read counts within the consensus regions #########################
replicates='Rep1 Rep2 Rep3 Rep4 Rep5'

#SH-SY5Y
for rep in $replicates
do
    featureCounts -a $indir/SY5Y_all_cat_btmerge.gff -t region -s 0 -g region_id -B -p -o $outdir/counts_SY5Y_${rep}_control_consensus.txt $indir/SY5Y_${rep}_control_ChIP.sorted.bam -T 10
done
wait

for rep in $replicates
do
    featureCounts -a $indir/SY5Y_all_cat_btmerge.gff -t region -s 0 -g region_id -B -p -o $outdir/counts_SY5Y_${rep}_5dPB_consensus.txt $indir/SY5Y_${rep}_5dPB_ChIP.sorted.bam -T 10
done
wait

#IMR-32
for rep in $replicates
do
    featureCounts -a $indir/IMR32_all_cat_btmerge.gff -t region -s 0 -g region_id -B -p -o $outdir/counts_IMR32_${rep}_control_consensus.txt $indir/IMR32_${rep}_control_ChIP.sorted.bam -T 10
done
wait

for rep in $replicates
do
    featureCounts -a $indir/IMR32_all_cat_btmerge.gff -t region -s 0 -g region_id -B -p -o $outdir/counts_IMR32_${rep}_5dPB_consensus.txt $indir/IMR32_${rep}_5dPB_ChIP.sorted.bam -T 12
done
wait

#SK-N-BE(2)C
for rep in $replicates
do
    featureCounts -a $indir/BE2C_all_cat_btmerge.gff -t region -s 0 -g region_id -B -p -o $outdir/counts_BE2C_${rep}_control_consensus.txt $indir/BE2C_${rep}_control_ChIP.sorted.bam -T 12
done
wait

for rep in $replicates
do
    featureCounts -a $indir/BE2C_all_cat_btmerge.gff -t region -s 0 -g region_id -B -p -o $outdir/counts_BE2C_${rep}_7dPB_consensus.txt $indir/BE2C_${rep}_7dPB_ChIP.sorted.bam -T 12
done
wait


############################## get library sizes ##############################

#SH-SY5Y
for rep in $replicates
do
    samtools flagstat $indir/SY5Y_${rep}_control_ChIP.sorted.bam > $outdir/SY5Y_${rep}_control_ChIP_flagstat.txt &
done
wait

for rep in $replicates
do
    samtools flagstat $indir/SY5Y_${rep}_5dPB_ChIP.sorted.bam > $outdir/SY5Y_${rep}_5dPB_ChIP_flagstat.txt &
done
wait

#IMR-32
for rep in $replicates
do
    samtools flagstat $indir/IMR32_${rep}_control_ChIP.sorted.bam > $outdir/IMR32_${rep}_control_ChIP_flagstat.txt &
done
wait

for rep in $replicates
do
    samtools flagstat $indir/IMR32_${rep}_5dPB_ChIP.sorted.bam > $outdir/IMR32_${rep}_5dPB_ChIP_flagstat.txt &
done
wait

#SK-N-BE(2)C
for rep in $replicates
do
    samtools flagstat $indir/BE2C_${rep}_control_ChIP.sorted.bam > $outdir/BE2C_${rep}_control_ChIP_flagstat.txt &
done
wait

for rep in $replicates
do
    samtools flagstat $indir/BE2C_${rep}_7dPB_ChIP.sorted.bam > $outdir/BE2C_${rep}_7dPB_ChIP_flagstat.txt &
done
wait


############################## adjust counts to normalise for library size ##############################
#SH-SY5Y
for rep in $replicates
do
    python sort_normalisation.py $outdir/counts_SY5Y_${rep}_control_consensus.txt $outdir/SY5Y_${rep}_control_ChIP_flagstat.txt $outdir/counts_SY5Y_${rep}_control_consensus_ENHANCER_REGION_MAP.txt SY5Y_${rep}_control &
done
wait

for rep in $replicates
do
    python sort_normalisation.py $outdir/counts_SY5Y_${rep}_5dPB_consensus.txt $outdir/SY5Y_${rep}_5dPB_ChIP_flagstat.txt $outdir/counts_SY5Y_${rep}_5dPB_consensus_ENHANCER_REGION_MAP.txt SY5Y_${rep}_5dPB &
done
wait

#IMR-32
for rep in $replicates
do
    python sort_normalisation.py $outdir/counts_IMR32_${rep}_control_consensus.txt $outdir/IMR32_${rep}_control_ChIP_flagstat.txt $outdir/counts_IMR32_${rep}_control_consensus_ENHANCER_REGION_MAP.txt IMR32_${rep}_control &
done
wait

for rep in $replicates
do
    python sort_normalisation.py $outdir/counts_IMR32_${rep}_5dPB_consensus.txt $outdir/IMR32_${rep}_5dPB_ChIP_flagstat.txt $outdir/counts_IMR32_${rep}_5dPB_consensus_ENHANCER_REGION_MAP.txt IMR32_${rep}_5dPB &
done
wait

#SK-N-BE(2)C
for rep in $replicates
do
    python sort_normalisation.py $outdir/counts_BE2C_${rep}_control_consensus.txt $outdir/BE2C_${rep}_control_ChIP_flagstat.txt $outdir/counts_BE2C_${rep}_control_consensus_ENHANCER_REGION_MAP.txt BE2C_${rep}_control &
done
wait

for rep in $replicates
do
    python sort_normalisation.py $outdir/counts_BE2C_${rep}_7dPB_consensus.txt $outdir/BE2C_${rep}_7dPB_ChIP_flagstat.txt $outdir/counts_BE2C_${rep}_7dPB_consensus_ENHANCER_REGION_MAP.txt BE2C_${rep}_7dPB &
done
wait

