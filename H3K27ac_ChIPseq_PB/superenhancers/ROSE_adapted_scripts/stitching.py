#!/usr/bin/env/ python
#stitching.py
#this script has been adapted from ROSE - doi: 10.1016/j.cell.2013.03.035 & doi: 10.1016/j.cell.2013.03.036

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
#=====================REGION STITCHING=============================
#==================================================================

def regionStitching(inputGFF,stitchWindow,tssWindow,annotFile,removeTSS):
    print('PERFORMING REGION STITCHING')
    #first have to turn bound region file into a locus collection
    #need to make sure this names correctly... each region should have a unique name
    boundCollection = ROSE_utils.gffToLocusCollection(inputGFF)

    debugOutput = []
    #filter out all bound regions that overlap the TSS of an ACTIVE GENE
    if removeTSS:
        #first make a locus collection of TSS
        startDict = ROSE_utils.makeStartDict(annotFile)

        #now makeTSS loci for active genes
        removeTicker=0
        #this loop makes a locus centered around +/- tssWindow of transcribed genes
        #then adds it to the list tssLoci
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,tssWindow,tssWindow))

        #this turns the tssLoci list into a LocusCollection
        #50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = ROSE_utils.LocusCollection(tssLoci,50)

        #gives all the loci in boundCollection
        boundLoci = boundCollection.getLoci()

        #this loop will check if each bound region is contained by the TSS exclusion zone
        #this will drop out a lot of the promoter only regions that are tiny
        #typical exclusion window is around 2kb
        for locus in boundLoci:
            if len(tssCollection.getContainers(locus,'both'))>0:
                
                #if true, the bound locus overlaps an active gene
                boundCollection.remove(locus)
                debugOutput.append([locus.__str__(),locus.ID(),'CONTAINED'])
                removeTicker+=1
        print('REMOVED %s LOCI BECAUSE THEY WERE CONTAINED BY A TSS' % (removeTicker))

    #boundCollection is now all enriched region loci that don't overlap an active TSS
    stitchedCollection = boundCollection.stitchCollection(stitchWindow,'both')

    if removeTSS:
        #now replace any stitched region that overlap 2 distinct genes
        #with the original loci that were there
        fixedLoci = []
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,50,50))


        #this turns the tssLoci list into a LocusCollection
        #50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = ROSE_utils.LocusCollection(tssLoci,50)
        removeTicker = 0
        originalTicker = 0
        for stitchedLocus in stitchedCollection.getLoci():
            overlappingTSSLoci = tssCollection.getOverlap(stitchedLocus,'both')
            tssNames = [startDict[tssLocus.ID()]['name'] for tssLocus in overlappingTSSLoci]
            tssNames = ROSE_utils.uniquify(tssNames)
            if len(tssNames) > 2:
            
                #stitchedCollection.remove(stitchedLocus)
                originalLoci = boundCollection.getOverlap(stitchedLocus,'both')
                originalTicker+=len(originalLoci)
                fixedLoci+=originalLoci
                debugOutput.append([stitchedLocus.__str__(),stitchedLocus.ID(),'MULTIPLE_TSS'])
                removeTicker+=1
            else:
                fixedLoci.append(stitchedLocus)

        print('REMOVED %s STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs' % (removeTicker))
        print('ADDED BACK %s ORIGINAL LOCI' % (originalTicker))
        fixedCollection = ROSE_utils.LocusCollection(fixedLoci,50)
        return fixedCollection,debugOutput
    else:
        return stitchedCollection,debugOutput


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
    parser.add_option("-t","--tss", dest="tss",nargs = 1, default=0,
                      help = "Enter a distance from TSS to exclude. 0 = no TSS exclusion")

    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()

    if not options.input or not options.out or not options.genome:
        print('missing inputs')
        parser.print_help()
        exit()

    #making the out folder if it doesn't exist
    outFolder = ROSE_utils.formatFolder(options.out,True)
    
    #figuring out folder schema
    gffFolder = ROSE_utils.formatFolder(outFolder+'gff/',True)
    mappedFolder = ROSE_utils.formatFolder(outFolder+ 'mappedGFF/',True)

    #GETTING INPUT FILE -> this should be gff anyway....
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
        print('WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT') #should probably just make this quit here
        #COPY THE INPUT GFF TO THE GFF FOLDER
        inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile,gffFolder))        

    #Stitch parameter
    stitchWindow = int(options.stitch)
    
    #tss options
    tssWindow = int(options.tss)
    if tssWindow != 0:
        removeTSS = True
        print('TSSexclusion is %s' % tssWindow)
    else:
        removeTSS = False
        print('TSSexclusion is %s' % tssWindow)

    #GETTING THE BOUND REGION FILE USED TO DEFINE ENHANCERS
    print('USING %s AS THE INPUT GFF' % (inputGFFFile))
    inputName = inputGFFFile.split('/')[-1].split('.')[0]

    #GETTING THE GENOME
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

    annotFile = genomeDict[upper(genome)]

    #MAKING THE START DICT
    print('MAKING START DICT')
    startDict = ROSE_utils.makeStartDict(annotFile)

    #LOADING IN THE BOUND REGION REFERENCE COLLECTION
    print('LOADING IN GFF REGIONS')
    referenceCollection = ROSE_utils.gffToLocusCollection(inputGFFFile)

    #NOW STITCH REGIONS
    print('STITCHING REGIONS TOGETHER')
    stitchedCollection,debugOutput = regionStitching(inputGFFFile,stitchWindow,tssWindow,annotFile,removeTSS)

    #NOW MAKE A STITCHED COLLECTION GFF
    print('MAKING GFF FROM STITCHED COLLECTION')
    stitchedGFF=ROSE_utils.locusCollectionToGFF(stitchedCollection)
    
    if not removeTSS:
        stitchedGFFFile = '%s%s_%sKB_STITCHED.gff' % (gffFolder,inputName,stitchWindow/1000)
        stitchedGFFName = '%s_%sKB_STITCHED' % (inputName,stitchWindow/1000)
        debugOutFile = '%s%s_%sKB_STITCHED.debug' % (gffFolder,inputName,stitchWindow/1000)
    else:
        stitchedGFFFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.gff' % (gffFolder,inputName,stitchWindow/1000)
        stitchedGFFName = '%s_%sKB_STITCHED_TSS_DISTAL' % (inputName,stitchWindow/1000)
        debugOutFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.debug' % (gffFolder,inputName,stitchWindow/1000)

    #WRITING DEBUG OUTPUT TO DISK
        
    if debug:
        print('WRITING DEBUG OUTPUT TO DISK AS %s' % (debugOutFile))
        ROSE_utils.unParseTable(debugOutput,debugOutFile,'\t')

    #WRITE THE GFF TO DISK
    print('WRITING STITCHED GFF TO DISK AS %s' % (stitchedGFFFile))
    ROSE_utils.unParseTable(stitchedGFF,stitchedGFFFile,'\t')



if __name__ == "__main__":
    main()
