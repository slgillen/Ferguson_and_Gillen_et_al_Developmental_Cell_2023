#!/usr/bin/env python 

import sys
fcounts = sys.argv[1]
fstat = sys.argv[2]
outname = sys.argv[3]
samplename = sys.argv[4]


import os
import subprocess
import math
import itertools


with open(fstat) as statfile:
    lines = statfile.readlines()
    for line in lines:
        if 'mapped (' in line:
            mapped_count=line.split(' ')[0]

statfile.close()

print(mapped_count)
mapped_count=math.floor(float(mapped_count)/2)
print(mapped_count)
scaling_factor=mapped_count/1000000
print(scaling_factor)


with open(fcounts) as infile:
    with open(outname, 'w') as outfile:
        header_names=['REGION_ID','CHROM','START','STOP','NUM_LOCI','CONSTITUENT_SIZE',samplename]
        outfile.write('\t'.join(header_names)+'\n')
        for line in itertools.islice(infile, 2, None):
            cols = line.split("\t")
            region_name=cols[0]
            num_loci=region_name.split('_')[0]
            raw_count=cols[6]
            norm_count=float(raw_count)/scaling_factor
            indices=[0,1,2,3]
            newcols=[cols[index] for index in indices]
            newline='\t'.join(newcols)+'\t'+str(num_loci)+'\t'+cols[5]+'\t'+str(norm_count)+'\n'
            outfile.write(newline)

infile.close()
outfile.close()

