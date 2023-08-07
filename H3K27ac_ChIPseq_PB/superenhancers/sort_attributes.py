#!/usr/bin/env python 

import sys
currentBED = sys.argv[1]
newGFF = sys.argv[2]


import os
import subprocess


with open(currentBED) as infile:
    with open(newGFF, 'w') as outfile:
        for line in infile:
            cols = line.split("\t")
            indices=[0,1,2]
            newcols=[cols[index] for index in indices]
            newline=cols[0]+'\t'+'ROSE_STITCHED_LOCI'+'\t'+'region'+'\t'
            newline=newline+cols[1]+'\t'+str(cols[2].strip('\n'))+'\t'+'.'+'\t'+'.'+'\t'+'.'+'\t'
            newcols[2]=str(newcols[2].strip('\n'))
            attr='_'.join(newcols)
            newline=newline+'region_id '+'\"'+attr+'\"'+'\n'
            outfile.write(newline)

infile.close()
outfile.close()

