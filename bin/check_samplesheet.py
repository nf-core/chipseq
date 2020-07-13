#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on April 4th 2019 to check nf-core/chipseq design file
#######################################################################
#######################################################################

import os
import sys
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Reformat nf-core/chipseq samplesheet and check its contents.'
Epilog = "Example usage: python check_samplesheet.py <DESIGN_FILE> <READ_MAPPING_FILE> <CONTROL_MAPPING_FILE>"

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('DESIGN_FILE', help="Input samplesheet/design file.")
argParser.add_argument('READ_MAPPING_FILE', help="Output design file containing sample ids and reads.")
argParser.add_argument('CONTROL_MAPPING_FILE', help="Output design file containing ip vs control mappings.")
args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def reformat_design(DesignFile,ReadMappingFile,ControlMappingFile):

    ERROR_STR = 'ERROR: Please check design file'
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2', 'antibody', 'control']

    ## CHECK HEADER
    fin = open(DesignFile,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
        sys.exit(1)

    numColList = []
    sampleMappingDict = {}
    antibodyDict = {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',')]
            group,replicate,fastQFiles,antibody,control = lspl[0],lspl[1],[x for x in lspl[2:-2] if x],lspl[-2],lspl[-1]

            ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
            numCols = len(lspl)
            if numCols not in [6]:
                print("{}: Invalid number of columns (should be 6)!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            numColList.append(numCols)

            ## CHECK GROUP ID DOESNT CONTAIN SPACES
            if group.find(' ') != -1:
                print("{}: Group id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            ## CHECK REPLICATE COLUMN IS INTEGER
            if not replicate.isdigit():
                print("{}: Replicate id not an integer!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            replicate = int(replicate)

            for fastq in fastQFiles:
                ## CHECK FASTQ FILE EXTENSION
                if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                    print("{}: FastQ file has incorrect extension (has to be '.fastq.gz' or 'fq.gz') - {}\nLine: '{}'".format(ERROR_STR,fastq,line.strip()))
                    sys.exit(1)

            ## CREATE GROUP MAPPING DICT = {GROUP_ID: {REPLICATE_ID:[[FASTQ_FILES]]}
            if group not in sampleMappingDict:
                sampleMappingDict[group] = {}
            if replicate not in sampleMappingDict[group]:
                sampleMappingDict[group][replicate] = []
            sampleMappingDict[group][replicate].append(fastQFiles)

            ## CHECK BOTH ANTIBODY AND CONTROL COLUMNS HAVE VALID VALUES
            if antibody:
                if antibody.find(' ') != -1:
                    print("{}: Antibody id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)
                if not control:
                    print("{}: both Antibody and Control must be specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)
            if control:
                if control.find(' ') != -1:
                    print("{}: Control id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)
                if not antibody:
                    print("{}: both Antibody and Control must be specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)

            ## CREATE ANTIBODY MAPPING CONTROL DICT
            if antibody and control:
                antibodyDict[group] = (antibody,control)

        else:
            fin.close()
            break

    ## CHECK IF DATA IS PAIRED-END OR SINGLE-END AND NOT A MIXTURE
    if min(numColList) != max(numColList):
        print("{}: Mixture of paired-end and single-end reads!".format(ERROR_STR))
        sys.exit(1)

    ## CHECK IF ANTIBODY AND CONTROL COLUMNS HAVE BEEN SPECIFIED AT LEAST ONCE
    if len(antibodyDict) == 0:
        print("{}: Antibody and Control must be specified at least once!".format(ERROR_STR))
        sys.exit(1)

    ## WRITE READ MAPPING FILE
    antibodyGroupDict = {}
    fout = open(ReadMappingFile,'w')
    fout.write(','.join(['sample_id','fastq_1','fastq_2']) + '\n')
    for group in sorted(sampleMappingDict.keys()):

        ## CHECK THAT REPLICATE IDS ARE IN FORMAT 1..<NUM_REPLICATES>
        uniq_rep_ids = set(sampleMappingDict[group].keys())
        if len(uniq_rep_ids) != max(uniq_rep_ids):
            print("{}: Replicate IDs must start with 1..<num_replicates>\nGroup: {}, Replicate IDs: {}".format(ERROR_STR,group,list(uniq_rep_ids)))
            sys.exit(1)

        ## RECONSTRUCT LINE FOR SAMPLE IN DESIGN
        for replicate in sorted(sampleMappingDict[group].keys()):
            for idx in range(len(sampleMappingDict[group][replicate])):
                fastQFiles = sampleMappingDict[group][replicate][idx]

                ## GET SAMPLE_ID,FASTQ_1,FASTQ_2 COLUMNS
                sample_id = "{}_R{}_T{}".format(group,replicate,idx+1)
                oList = [sample_id] + fastQFiles
                if len(fastQFiles) == 1:
                    oList += ['']
                fout.write(','.join(oList) + '\n')

                ## EXTRAPOLATE CONTROL COLUMN
                if group in antibodyDict:
                    antibody,control = antibodyDict[group]
                    if control in sampleMappingDict.keys():
                        control_id = "{}_R1".format(control)
                        if replicate in sampleMappingDict[control]:
                            control_id = "{}_R{}".format(control,replicate)
                        if antibody not in antibodyGroupDict:
                            antibodyGroupDict[antibody] = {}
                        if group not in antibodyGroupDict[antibody]:
                            antibodyGroupDict[antibody][group] = []
                        antibodyList = [sample_id[:-3],control_id]
                        if not antibodyList in antibodyGroupDict[antibody][group]:
                            antibodyGroupDict[antibody][group].append(antibodyList)
                    else:
                        print("{}: Control id not a valid group\nControl id: {}, Valid Groups: {}".format(ERROR_STR,control,sorted(sampleMappingDict.keys())))
                        sys.exit(1)
    fout.close()

    ## WRITE SAMPLE TO CONTROL MAPPING FILE
    fout = open(ControlMappingFile,'w')
    fout.write(','.join(['sample_id','control_id','antibody','replicatesExist','multipleGroups']) + '\n')
    for antibody in sorted(antibodyGroupDict.keys()):
        repsExist = '0'
        if max([len(x) for x in antibodyGroupDict[antibody].values()]) > 1:
            repsExist = '1'
        multipleGroups = '0'
        if len(antibodyGroupDict[antibody].keys()) > 1:
            multipleGroups = '1'
        for group in sorted(antibodyGroupDict[antibody].keys()):
            for antibodyList in antibodyGroupDict[antibody][group]:
                fout.write(','.join(antibodyList+[antibody,repsExist,multipleGroups]) + '\n')
    fout.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

reformat_design(DesignFile=args.DESIGN_FILE,ReadMappingFile=args.READ_MAPPING_FILE,ControlMappingFile=args.CONTROL_MAPPING_FILE)

############################################
############################################
############################################
############################################
