#!/bin/sh

set -e

# IN-FILES
if [ -z "$1" ]
then
    echo usage $0 \<alignment\> \<annotation_gff\> \<comparison_outfile\> \<coverage_outfile\>
    exit 2
else
    alignment=$1
fi

if [ -z "$2" ]
then
    echo usage $0 \<alignment\> \<annotation_gff\> \<comparison_outfile\> \<coverage_outfile\>
    exit 2
else
    annotation=$2
fi

# OUT-FILES
if [ -z "$3" ]
then
    echo usage $0 \<alignment\> \<annotation_gff\> \<comparison_outfile\> \<coverage_outfile\>
    exit 2
else
    comparison_out=$3
fi

if [ -z "$4" ]
then
    echo usage $0 \<alignment\> \<annotation_gff\> \<comparison_outfile\> \<coverage_outfile\>
    exit 2
else
    coverage_out=$4
fi

# INTERMEDIATE-FILES
tmpoutdir=`mktemp -d` || exit 1
annotation_intron_list=$tmpoutdir/annotation.introns

echo "Create intron list from annotation:"
echo "python ./tools/gen_intronlist_from_annotation.py -a $annotation -o $annotation_intron_list"
python ./tools/gen_intronlist_from_annotation.py -a $annotation -o $annotation_intron_list
echo ""

echo "Compare spliced alignment to annotated introns:"
echo "python compare_intron_lists.py -O 20 -I $annotation_intron_list -C ${annotation_intron_list}.cov -a $alignment -H $comparison_out -g $coverage_out"
python compare_intron_lists.py -O 20 -I $annotation_intron_list -C ${annotation_intron_list}.cov -a $alignment -H $comparison_out -g $coverage_out

# Clean up
rm -r $tmpoutdir
