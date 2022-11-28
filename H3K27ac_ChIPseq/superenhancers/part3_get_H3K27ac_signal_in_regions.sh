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
