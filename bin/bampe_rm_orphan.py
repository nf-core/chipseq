#!/usr/bin/env python

###############################################################################
###############################################################################
## Created on February 1st 2017 to remove singletons from paired-end BAM file
###############################################################################
###############################################################################

import os
import pysam
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Remove singleton reads from paired-end BAM file i.e if read1 is present in BAM file without read 2 and vice versa.'
Epilog = """Example usage: bampe_rm_orphan.py <BAM_INPUT_FILE> <BAM_OUTPUT_FILE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('BAM_INPUT_FILE', help="Input BAM file sorted by name.")
argParser.add_argument('BAM_OUTPUT_FILE', help="Output BAM file sorted by name.")

## OPTIONAL PARAMETERS
argParser.add_argument('-fr', '--only_fr_pairs', dest="ONLY_FR_PAIRS", help="Only keeps pairs that are in FR orientation on same chromosome.",action='store_true')
args = argParser.parse_args()

############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def makedir(path):

    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def bampe_rm_orphan(BAMIn,BAMOut,onlyFRPairs=False):

    ## SETUP DIRECTORY/FILE STRUCTURE
    OutDir = os.path.dirname(BAMOut)
    makedir(OutDir)

    ## COUNT VARIABLES
    totalReads = 0; totalOutputPairs = 0; totalSingletons = 0; totalImproperPairs = 0

    ## ITERATE THROUGH BAM FILE
    EOF = 0
    SAMFin = pysam.Samfile(BAMIn,"rb")
    SAMFout = pysam.Samfile(BAMOut, "wb",header=SAMFin.header)
    currRead = SAMFin.next()
    for read in SAMFin:
        totalReads += 1
        if currRead.qname == read.qname:
            pair1 = currRead; pair2 = read

            ## FILTER FOR READS ON SAME CHROMOSOME IN FR ORIENTATION
            if onlyFRPairs:
                if pair1.tid == pair2.tid:

                    ## READ1 FORWARD AND READ2 REVERSE STRAND
                    if not pair1.is_reverse and pair2.is_reverse:
                        if pair1.reference_start <= pair2.reference_start:
                            totalOutputPairs += 1
                            SAMFout.write(pair1)
                            SAMFout.write(pair2)
                        else:
                            totalImproperPairs += 1

                    ## READ1 REVERSE AND READ2 FORWARD STRAND
                    elif pair1.is_reverse and not pair2.is_reverse:
                        if pair2.reference_start <= pair1.reference_start:
                            totalOutputPairs += 1
                            SAMFout.write(pair1)
                            SAMFout.write(pair2)
                        else:
                            totalImproperPairs += 1

                    else:
                        totalImproperPairs += 1
                else:
                    totalImproperPairs += 1
            else:
                totalOutputPairs += 1
                SAMFout.write(pair1)
                SAMFout.write(pair2)

            ## RESET COUNTER
            try:
                totalReads += 1
                currRead = SAMFin.next()
            except:
                StopIteration
                EOF = 1

        ## READS WHERE ONLY ONE OF A PAIR IS IN FILE
        else:
            totalSingletons += 1
            pair1 = currRead
            currRead = read

    if not EOF:
        totalReads += 1
        totalSingletons += 1
        pair1 = currRead

    ## CLOSE ALL FILE HANDLES
    SAMFin.close()
    SAMFout.close()

    LogFile = os.path.join(OutDir,'%s_bampe_rm_orphan.log' % (os.path.basename(BAMOut[:-4])))
    SamLogFile = open(LogFile,'w')
    SamLogFile.write('\n##############################\n')
    SamLogFile.write('FILES/DIRECTORIES')
    SamLogFile.write('\n##############################\n\n')
    SamLogFile.write('Input File: ' + BAMIn + '\n')
    SamLogFile.write('Output File: ' + BAMOut + '\n')
    SamLogFile.write('\n##############################\n')
    SamLogFile.write('OVERALL COUNTS')
    SamLogFile.write('\n##############################\n\n')
    SamLogFile.write('Total Input Reads = ' + str(totalReads) + '\n')
    SamLogFile.write('Total Output Pairs = ' + str(totalOutputPairs) + '\n')
    SamLogFile.write('Total Singletons Excluded = ' + str(totalSingletons) + '\n')
    SamLogFile.write('Total Improper Pairs Excluded = ' + str(totalImproperPairs) + '\n')
    SamLogFile.write('\n##############################\n')
    SamLogFile.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

bampe_rm_orphan(BAMIn=args.BAM_INPUT_FILE,BAMOut=args.BAM_OUTPUT_FILE,onlyFRPairs=args.ONLY_FR_PAIRS)

############################################
############################################
############################################
############################################
