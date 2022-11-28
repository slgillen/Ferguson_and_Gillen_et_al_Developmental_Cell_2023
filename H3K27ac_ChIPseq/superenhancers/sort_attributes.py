#!/usr/bin/env python 

import sys
currentGFF = sys.argv[1]
newGFF = sys.argv[2]


import os
import subprocess


with open(currentGFF) as infile:
    with open(newGFF, 'w') as outfile:
        for line in infile:
            cols = line.split("\t")
            #NOTE: there are empty columns
            indices=[3,4]
            newcols=[cols[index] for index in indices]
            newline=cols[0]+'\t'+'ROSE_STITCHED_LOCI'+'\t'+'region'+'\t'
            newline=newline+'\t'.join(newcols)+'\t'+'.'+'\t'+'.'+'\t'+'.'+'\t'+'region_id '+'\"'+str(cols[8].strip('\n'))+'\"'+'\n'
            outfile.write(newline)

infile.close()
outfile.close()

