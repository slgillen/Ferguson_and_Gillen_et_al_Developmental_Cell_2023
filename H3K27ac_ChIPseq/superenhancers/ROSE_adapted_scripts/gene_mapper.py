#!/bin/bash
#adapted from ROSE
'''
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
APRIL 11, 2013
VERSION 0.1
CONTACT: youngcomputation@wi.mit.edu
'''

import sys
import ROSE_utils
import time
import os
from string import upper,join
from collections import defaultdict

#==================================================================
#=========================MAIN METHOD==============================
#==================================================================

def main():
    '''
    main run call
    '''
    debug = False


    from optparse import OptionParser
    usage = "usage: %prog [options] -g [GENOME] -i [INPUT_REGION_GFF] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-i","--i", dest="input",nargs = 1, default=None,
                      help = "Enter a .gff or .bed file of binding sites used to make enhancers")
    parser.add_option("-o","--out", dest="out",nargs = 1, default=None,
                      help = "Enter an output folder")
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the genome build (MM9,MM8,HG18,HG19)")
    
    #optional flags
    parser.add_option("-s","--stitch", dest="stitch",nargs = 1, default=12500,
                      help = "Enter a max linking distance for stitching")

    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()


    if not options.input or not options.out or not options.genome:
        print('missing input')
        parser.print_help()
        exit()

    #making the out folder if it doesn't exist
    outFolder = ROSE_utils.formatFolder(options.out,True)

    #figuring out folder schema
    gffFolder = ROSE_utils.formatFolder(outFolder+'gff/',True)
    mappedFolder = ROSE_utils.formatFolder(outFolder+ 'mappedGFF/',True)


    #GETTING INPUT FILE
    if options.input.split('.')[-1] == 'bed':
        #CONVERTING A BED TO GFF
        inputGFFName = options.input.split('/')[-1][0:-4]
        inputGFFFile = '%s%s.gff' % (gffFolder,inputGFFName)
        ROSE_utils.bedToGFF(options.input,inputGFFFile)
    elif options.input.split('.')[-1] =='gff':
        #COPY THE INPUT GFF TO THE GFF FOLDER
        inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile,gffFolder))        

    else:
        print('WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT')
        #COPY THE INPUT GFF TO THE GFF FOLDER
        inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile,gffFolder))        

    #GETTING THE LIST OF BAMFILES TO PROCESS
    #if options.bams:
    #    bamFileList = options.bams.split(',')
    #    bamFileList = ROSE_utils.uniquify(bamFileList)
    
    

    #GETTING THE BOUND REGION FILE USED TO DEFINE ENHANCERS
    print('USING %s AS THE INPUT GFF' % (inputGFFFile))
    inputName = inputGFFFile.split('/')[-1].split('.')[0]

    #GETING THE GENOME
    genome = options.genome
    print('USING %s AS THE GENOME' % genome)

    #GETTING THE CORRECT ANNOT FILE
    cwd = os.getcwd()
    genomeDict = {
        'HG18':'%s/annotation/hg18_refseq.ucsc' % (cwd),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (cwd),
        'HG19':'%s/annotation/hg19_refseq.ucsc' % (cwd),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (cwd),
        'MM10':'%s/annotation/mm10_refseq.ucsc' % (cwd),
        }

    #superTableFile = "%s_SuperEnhancers.table.txt" % ('calledSupers/'+inputName)
    #superTableFile = "%s.bed" % ('superDepthOverlap/'+inputName)
    #superTableFile = "%s.bed" % ('manually_calledSupers/'+inputName)
    superTableFile = "%s.bed" % (inputName)
    print('OUTFOLDER is %s' % outFolder)
    cmd = "python ROSE_geneMapper.py -o %s -g %s -i %s%s" % (outFolder+'mapped_genes/',genome,outFolder,superTableFile)
    os.system(cmd)    




if __name__ == "__main__":
    main()
