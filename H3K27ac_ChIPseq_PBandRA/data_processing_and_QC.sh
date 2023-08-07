#!/bin/bash

########  tools used in this script ########
#fastp version: 0.12.4
#fastqc version: 0.11.9
#bowtie2 version: 2.4.4
#bedtools version: 2.26.0
#samtools version: 1.9
#MACS2 version: 2.2.7.1
#deeptools version: 3.5.1
#sambamba version: 0.6.6

###### directories ######
maindir='' #main directory of inputs and outputs
indir=$maindir/raw_fastq
indexdir=$maindir/indexes

#create output directories
fastqcdir=$maindir/fastqc
mkdir $fastqcdir
trimdir=$maindir/fastp
mkdir $trimdir
fastqcdir2=$maindir/fastqc_after_fastp
mkdir $fastqcdir2
aligndir=$maindir/bowtie2_aligned
mkdir $aligndir
filtdir=$maindir/bam_filtering
mkdir $filtdir
bwdir=$maindir/bigwig
mkdir $bwdir


####################### read in all samples names ##########################
#assumes run from directory containing the fastq files
samplenames=$(find . -name BE2C\*\_R1\.fastq\.gz | sed 's/.\///' | sed s/\_R1\.fastq\.gz//)
echo $samplenames


####################### initial Q.C check ##########################
for s in ${samplenames[@]}
do
  fastqc $indir/${s}_R1.fastq.gz --outdir=$fastqcdir &
done
wait

for s in ${samplenames[@]}
do
  fastqc $indir/${s}_R2.fastq.gz --outdir=$fastqcdir &
done
wait

echo "initial fastqc complete"

####################### data trimming with fastp ##########################
for s in ${samplenames[@]}
do	
  fastp -i $indir/${s}_R1.fastq.gz -I $indir/${s}_R2.fastq.gz --trim_poly_g -w 8 -q 20 -l 20 -o $trimdir/${s}_fastp_R1.fq.gz -O $trimdir/${s}_fastp_R2.fq.gz > $trimdir/${s}_fastp_log.txt &
done
wait

echo "fastp complete"

####################### Q.C check after fastp trimming ##########################
for s in ${samplenames[@]}
do
  fastqc $trimdir/${s}_fastp_R1.fq.gz --outdir=$fastqcdir2 &
done
wait

for s in ${samplenames[@]}
do
  fastqc $trimdir/${s}_fastp_R2.fq.gz --outdir=$fastqcdir2 &
done
wait

echo "fastqc after fastp complete"


####################### alignment to hg19 with bowtie2 ##########################
#indexes made with bowtie2-build --threads 25 hg19.fa hg19

for s in ${samplenames[@]}
do
    (bowtie2 -p 60 -q -N 0 --dovetail --fr --no-mixed --no-unal -x $indexdir/hg19 -1 $trimdir/${s}_fastp_R1.fq.gz -2 $trimdir/${s}_fastp_R2.fq.gz -X 1000 --fr -S $aligndir/${s}.sam) 2> $aligndir/bowtie2_stats_${s}.txt
done
wait

####################### convert sam to bam ##########################
for s in ${samplenames[@]}
do
    samtools view -@ 20 -bS $aligndir/${s}.sam > $aligndir/${s}.bam
done
wait

####################### remove sam ##########################
for s in ${samplenames[@]}
do
    rm $aligndir/${s}.sam
done
wait

####################### Q.C check for alignment stats etc ##########################
for s in ${samplenames[@]}
do
    samtools flagstat --threads 20 $aligndir/${s}.bam > $aligndir/flagstat_${s}.txt
done
wait

####################### bam filtering ##########################
for s in ${samplenames[@]}
do
    sambamba view -h -t 20 -f bam -F "[XS] == null and not unmapped and proper_pair and not duplicate" $aligndir/${s}.bam > $filtdir/${s}_filt.bam
done
wait

####################### sort and index the bam file and remove unsorted file ##########################
#sort
for s in ${samplenames[@]}
do	
  samtools sort -@ 20 -o $filtdir/${s}_filt_sorted.bam $filtdir/${s}_filt.bam
done
wait

#remove unsorted
for s in ${samplenames[@]}
do	
  rm $filtdir/${s}_filt.bam
done
wait

#index
for s in ${samplenames[@]}
do	
  samtools index $filtdir/${s}_filt_sorted.bam $filtdir/${s}_filt_sorted.bai &
done
wait


####################### check paired-end fragment size distribution ##########################
for s in ${samplenames[@]}
do
    bamPEFragmentSize -p 20 -b $filtdir/${s}_filt_sorted.bam -hist $filtdir/PEfragsizes_${s}.png > $filtdir/fragsizeinfo_${s}.txt
done
wait

####################### flagstat check on filtered bam files ##########################
for s in ${samplenames[@]}
do
    samtools flagstat --threads 20 $filtdir/${s}_filt_sorted.bam > $filtdir/flagstat_filtcheck_${s}.txt
done
wait


####################### convert bam to bigwig ##########################

for s in ${samplenames[@]}
do
  echo $file
  bamCoverage -b $filtdir/${s}_filt_sorted.bam -o $bwdir/${s}_CPM.bw --binSize 10 -p 6 --normalizeUsing CPM --effectiveGenomeSize 2864785220 --extendReads
done
wait


####################### peak calling ##########################
for s in ${samplenames[@]}
do	
  echo ${s}
  macs2 callpeak -t $filtdir/${s}_IP_filt_sorted.bam -c $filtdir/${s}_Input_filt_sorted.bam -g 2.86e9 -f BAMPE --keep-dup all --broad --outdir $peakdir -n ${s} -B 2> $peakdir/${s}_broad_macs2.log
done
wait




