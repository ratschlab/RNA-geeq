#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Andre Kahles
# Copyright (C) 2009-2011 Max Planck Society
#

set -e 

PROG=`basename $0`

echo
echo ${PROG}: This program is part of the SAFT version 0.2.
echo
echo SAFT determines optimal filtering parameters for a given 
echo alignment based on a given annotation
echo 

if [ -z "$1" -o "$1" == '--help' ];
then
  echo Usage: $0 small\|big
  echo "   or:" $0 --help
  false
fi 
if [ "$1" != 'small' -a "$1" != 'big' ];
then
  echo invalid parameter
  false
fi

if [ "$1" == 'small' ];
then
  echo Note: Running this script takes about 1 minute \(on a single CPU\).
  GFF3_INPUT=data/nGASP-Train-I.gff3
  SAM_INPUT=data/nGASP-Train-I.sam
  BAM_INPUT=data/nGASP-Train-I.bam
  EXP=nGASP-Train-I
fi
if [ "$1" == 'big' ];
then
  echo Note: Running this script takes about 5 minutes \(on a single CPU\).
  GFF3_INPUT=data/nGASP-Train.gff3
  SAM_INPUT=data/nGASP-Train.sam
  BAM_INPUT=data/nGASP-Train.bam
  EXP=nGASP-Train
fi
LOAD_PROFILES="0"
PROFILES_FN="dummy_pr_in"
LEARN_PROFILES="0"
PROFILES_FN_OUT="dummy_pr_out"

RESULTDIR=./results-$1
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

#echo convert the alignments in SAM format to BAM format
#../tools/./sam_to_bam.sh ${SAMTOOLS_DIR} data ${FASTA_INPUT} ${SAM_INPUT}
#
SAFT_RES_DIR=$RESULTDIR
mkdir -p $SAFT_RES_DIR

OUT_BASE="${SAFT_RES_DIR}/`basename $SAM_INPUT | sed -e \"s/.sam$//g\"`"

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Generating annotation intron list %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo "Generating annotation intron list"
echo ""

if [ ! -f ${SAFT_RES_DIR}/`basename $GFF3_INPUT`.introns ]
then
    python ../gen_intronlist_from_annotation.py -v -a $GFF3_INPUT -o ${SAFT_RES_DIR}/`basename $GFF3_INPUT`.introns
fi

echo "done"
echo ""


echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Generating alignment feature list %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo 

if [ ! -f ${OUT_BASE}.features ]
then
    python ../get_intron_features.py -v -a $SAM_INPUT -o ${OUT_BASE}.features
fi

echo "done"
echo ""

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Search for optimal filter setting %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo 
python ../find_optimal_param_set.py -v -i ${SAFT_RES_DIR}/`basename $GFF3_INPUT`.introns -f ${OUT_BASE}.features -b ${OUT_BASE}.best_score -m ${OUT_BASE}.scoring_matrix
echo "done"
echo ""

if [ ! -z "$filter_out" ]
then

    echo
    echo %%%%%%%%%%%%%%%%%%%%%%%
    echo % 4. Filter Alignment % 
    echo %%%%%%%%%%%%%%%%%%%%%%%
    echo 
    min_ex_len=`tail -n 1 ${OUT_BASE}.best_score | cut -f 1`
    max_mm=`tail -n 1 ${OUT_BASE}.best_score | cut -f 2`
    min_support=`tail -n 1 ${OUT_BASE}.best_score | cut -f 3`
    support_string=""

    if [ "$min_support" != "1" ]
    then
        python ../filter_features.py -e $min_ex_len -X $max_mm -m $min_support -i ${OUT_BASE}.features -o ${OUT_BASE}.features_filtered
        support_string="-i ${OUT_BASE}.features_filtered"
    fi

    python ../filter_alignment.py -a $SAM_INPUT -o ${OUT_BASE}.filtered.sam -e $min_ex_len -X $max_mm $support_string
    echo ""
fi

echo
echo Filtering result is now available in ${OUT_BASE}.filtered.sam
echo

echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
