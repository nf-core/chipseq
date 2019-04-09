#!/usr/bin/env bash

## ARGUMENTS
INPUT_DIR=$1	  ## Path to input directory
FILE_GROUP=$2 	## File group identified
FILE_PATH=$3    ## Path to results/ directory containing files

## CONSENSUS PEAKS
find -L $INPUT_DIR -type f -iname "consensus*$FILE_GROUP*" | cut -c 2- | awk -v DIR="$FILE_PATH" -v OFS='\t' '{ print DIR$1, "0,0,0" }'

## DIFFERENTIAL INTERVALS FOR FDR <= 0.01
for fpath in $(find -L $INPUT_DIR -type f -iname "*$FILE_GROUP*" | grep "FDR" | grep -Ev "FDR0.01" | sort | cut -c 2-)
    do
        fdir=$(dirname "${fpath}")
        fbase=$(basename "${fpath}")
	      odir=$(echo "${fbase}" | awk -v var="."$FILE_GROUP 'BEGIN {FS=var;}{print $1}')
        opath=$FILE_PATH/$fdir/$odir/$fbase
        awk -v var="$opath" -v OFS='\t' 'BEGIN {print var, "255,0,0"}'
    done

## PEAK AND BIGWIG FILES
find -L $INPUT_DIR -type f -iname "*$FILE_GROUP*" | grep -Ev "consensus|FDR" | sort | cut -c 2- | awk -v DIR="$FILE_PATH" -v OFS='\t' '{ print DIR$1, "0,0,178" }'
