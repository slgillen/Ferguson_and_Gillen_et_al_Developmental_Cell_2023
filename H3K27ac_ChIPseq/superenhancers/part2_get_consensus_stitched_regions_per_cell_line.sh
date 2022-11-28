#!/bin/bash

#assumes running from within directory containing the gffs of the stitched regions


######################## sort file format ########################
Rscript convert_file_format.R

################### get consensus stitched region set across replicates and conditions per cell line ###################

maindir='' #base directory
mergedir=$maindir/merged_stitched_regions
mkdir $mergedir

#concatenate all situations
cat BE2C*.bed > $mergedir/BE2C_all_cat.bed
cat IMR32*.bed > $mergedir/IMR32_all_cat.bed
cat SY5Y*.bed > $mergedir/SY5Y_all_cat.bed

#sort
sort -k1,1 -k2,2n $mergedir/BE2C_all_cat.bed > $mergedir/BE2C_all_cat.sorted.bed
sort -k1,1 -k2,2n $mergedir/IMR32_all_cat.bed > $mergedir/IMR32_all_cat.sorted.bed
sort -k1,1 -k2,2n $mergedir/SY5Y_all_cat.bed > $mergedir/SY5Y_all_cat.sorted.bed

#this is the merged potential super-enhancer regions per cell line
#bedtools merge #default -d 0 i.e. overlapping and bookended merging
bedtools merge -i $mergedir/BE2C_all_cat.sorted.bed > $mergedir/BE2C_all_cat_btmerge.bed
bedtools merge -i $mergedir/IMR32_all_cat.sorted.bed > $mergedir/IMR32_all_cat_btmerge.bed
bedtools merge -i $mergedir/SY5Y_all_cat.sorted.bed > $mergedir/SY5Y_all_cat_btmerge.bed

#get number of regions in each consensus
wc -l $mergedir/BE2C_all_cat_btmerge.bed 
wc -l $mergedir/IMR32_all_cat_btmerge.bed 
wc -l $mergedir/SY5Y_all_cat_btmerge.bed 
