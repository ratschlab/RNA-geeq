#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    annotation="$1"
fi

if [ -z "$2" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    alignment="$2"
fi

if [ -z "$3" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    best_score="$3"
fi

if [ -z "$4" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    score_matrix="$4"
fi

if [ -z "$5" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    filter_out="$5"
fi

echo "Generating annotation intron list"
echo ""

if [ ! -f ${annotation}.introns ]
then
    python gen_intronlist_from_annotation.py -v -a $annotation -o ${annotation}.introns
fi

echo "done"
echo ""


echo "Generating alignment feature list"
echo ""

if [ ! -f ${alignment}.features ]
then
    python get_intron_features.py -v -a $alignment -o ${alignment}.features
fi

echo "done"
echo ""

echo "Search for optimal filter setting"
echo ""
python find_optimal_param_set.py -v -i ${annotation}.introns -f ${alignment}.features -b ${alignment}.best_score -m ${alignment}.scoring_matrix
echo "done"
echo ""

if [ ! -z "$filter_out" ]
then
    echo "Filter Alignment"
    echo ""
    min_ex_len=`tail -n 1 ${alignment}.best_score | cut -f 1`
    max_mm=`tail -n 1 ${alignment}.best_score | cut -f 2`
    min_support=`tail -n 1 ${alignment}.best_score | cut -f 3`
    support_string=""

    if [ "$min_support" != "1" ]
    then
        python filter_features.py -e $min_ex_len -X $max_mm -m $min_support -i ${alignment}.features -o ${alignment}.features_filtered
        support_string="-i ${alignment}.features_filtered"
    fi

    python filter_alignement.py -b -a $alignment -o $filter_out -e $min_ex_len -X $max_mm $support_string -s `which samtools`
    echo ""
fi

