#!/bin/bash

#this script uses the ROSE gene mapper to assign genes within 100kb to the superenhancers

superdir=''
outdir=$maindir/mapped_genes 
mkdir $outdir

python gene_mapper.py -g HG19 -i $superdir/BE2C_supertable.bed -o $superdir/ 
python gene_mapper.py -g HG19 -i $superdir/IMR32_supertable.bed -o $superdir/ 
python gene_mapper.py -g HG19 -i $superdir/SY5Y_supertable.bed -o $superdir/ 
