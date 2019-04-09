#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on April 4th 2019 to reformat nf-core/chipseq design file
#######################################################################
#######################################################################

import os
import sys
import requests
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Reformat nf-core/chipseq design file and check its contents.'
Epilog = """Example usage: python reformat_design.py <DESIGN_FILE_IN> <DESIGN_FILE_OUT>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('DESIGN_FILE_IN', help="Input design file.")
argParser.add_argument('DESIGN_FILE_OUT', help="Output design file.")
args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def reformat_design(DesignFileIn,DesignFileOut):

    ERROR_STR = 'ERROR: Please check design file'
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2', 'control']

    ## CHECK HEADER
    fin = open(DesignFileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print "{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER))
        sys.exit(1)

    numColList = []
    groupRepDict = {}
    groupControlDict = {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',')]

            ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
            numCols = len(lspl)
            if numCols not in [5]:
                print "{}: Invalid number of columns (should be 5)!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)
            numColList.append(numCols)

            ## CHECK GROUP COLUMN HAS NO SPACES
            group,replicate,fastQFiles,control = lspl[0],lspl[1],[x for x in lspl[2:-1] if x],lspl[-1]
            if group.find(' ') != -1:
                print "{}: Group id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)

            ## CHECK REPLICATE COLUMN IS INTEGER
            if not replicate.isdigit():
                print "{}: Replicate id not an integer!\nLine: '{}'".format(ERROR_STR,line.strip())
                sys.exit(1)

            for fastq in fastQFiles:
                ## CHECK FASTQ FILE EXTENSION
                if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                    print "{}: FastQ file has incorrect extension (has to be '.fastq.gz' or 'fq.gz') - {}\nLine: '{}'".format(ERROR_STR,fastq,line.strip())
                    sys.exit(1)

                ## CHECK FASTQ FILES EXIST PER SAMPLE
                if fastq[:4] not in ['http']:
                    if not os.path.exists(fastq):
                        print "{}: FastQ file does not exist - {}\nLine: '{}'".format(ERROR_STR,fastq,line.strip())
                        sys.exit(1)
                else:
                    if requests.head(fastq).status_code >= 400:
                        print "{}: FastQ file does not exist - {}\nLine: '{}'".format(ERROR_STR,fastq,line.strip())
                        sys.exit(1)

            ## CREATE GROUP MAPPING DICT = {GROUP_ID: {REPLICATE_ID:[[FASTQ_FILES]]}
            replicate = int(replicate)
            if not groupRepDict.has_key(group):
                groupRepDict[group] = {}
            if not groupRepDict[group].has_key(replicate):
                groupRepDict[group][replicate] = []
            groupRepDict[group][replicate].append(fastQFiles)

            ## CREATE GROUP CONTROL DICT = {GROUP_ID: CONTROL_ID}
            if control:
                groupControlDict[group] = control

        else:
            fin.close()
            break

    ## CHECK IF DATA IS PAIRED-END OR SINGLE-END AND NOT A MIXTURE
    if min(numColList) != max(numColList):
        print "{}: Mixture of paired-end and single-end reads!".format(ERROR_STR)
        sys.exit(1)

    ## CHECK IF MULTIPLE GROUPS EXIST
    multiGroups = False
    if len(groupRepDict) > 1:
        multiGroups = True

    ## WRITE TO FILE
    numRepList = []
    fout = open(DesignFileOut,'w')
    fout.write(','.join(['sample_id','fastq_1','fastq_2','control_id']) + '\n')
    for group in sorted(groupRepDict.keys()):

        ## CHECK THAT REPLICATE IDS ARE IN FORMAT 1..<NUM_REPLICATES>
        uniq_rep_ids = set(groupRepDict[group].keys())
        if len(uniq_rep_ids) != max(uniq_rep_ids):
            print "{}: Replicate IDs must start with 1..<num_replicates>\nGroup: {}, Replicate IDs: {}".format(ERROR_STR,group,list(uniq_rep_ids))
            sys.exit(1)
        numRepList.append(max(uniq_rep_ids))

        ## RECONSTRUCT LINE FOR SAMPLE IN DESIGN
        for replicate in sorted(groupRepDict[group].keys()):
            for idx in range(len(groupRepDict[group][replicate])):
                fastQFiles = groupRepDict[group][replicate][idx]

                ## GET SAMPLE_ID,FASTQ_1,FASTQ_2 COLUMNS
                sample_id = "{}_R{}_T{}".format(group,replicate,idx+1)
                oList = [sample_id] + fastQFiles
                if len(fastQFiles) == 1:
                    oList += ['']

                ## EXTRAPOLATE CONTROL COLUMN
                control_col = ''
                if groupControlDict.has_key(group):
                    control = groupControlDict[group]
                    if control in groupRepDict.keys():
                        if groupRepDict[control].has_key(replicate):
                            control_col += "{}_R{}".format(control,replicate)
                        else:
                            control_col += "{}_R1".format(control)
                    else:
                        print "{}: Control id not a valid group\nControl id: {}, Valid Groups: {}".format(ERROR_STR,groupControlDict[group],sorted(groupRepDict.keys()))
                        sys.exit(1)
                oList += [control_col]
                fout.write(','.join(oList) + '\n')
    fout.close()

    ## CHECK IF REPLICATES IN DESIGN
    repsExist = False
    if max(numRepList) != 1:
        repsExist = True

    ## CHECK FOR BALANCED DESIGN ACROSS MULTIPLE GROUPS.
    balancedDesign = False
    if len(set(numRepList)) == 1 and multiGroups and repsExist:
        balancedDesign = True

############################################
############################################
## RUN FUNCTION
############################################
############################################

reformat_design(DesignFileIn=args.DESIGN_FILE_IN,DesignFileOut=args.DESIGN_FILE_OUT)

############################################
############################################
############################################
############################################
