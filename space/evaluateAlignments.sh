#!/bin/sh

set -e

# IN-FILES
alignment=$1
annotation=$4
genome_fasta=$2
# INTERMEDIATE-FILES
tmpoutdir=`mktemp -d` || exit 1
annotation_intron_list=$tmpoutdir/annotation.introns

# OUT-FILES
statistics_out=
comparison_out=
coverage_out=

if [ "$#" == "0" ];
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<statistics_outfile\> [\<annotation_gff\> \<comparison_outfile\> \<coverage_outfile\>]
        exit -1
fi

# Check for genome_config
if [ ! -d "$tmpoutdir/genome" ];
then
    echo "creating genome information object from $genome_fasta in $tmpoutdir ..."
        python make_gio.py $genome_fasta $tmpoutdir
    echo "... done"
    echo ""
fi

echo "Generate statistics:"
echo "python gen_alignment_statistics.py -a $alignment -g $genome_config -H $statistics_out"
python gen_alignment_statistics.py -a $alignment -g $genome_config -H $statistics_out
echo ""

# OPTIONAL FOR SPLICED ALIGNMENTS

echo "Create intron list from annotation:"
echo "python gen_intronlist_from_annotation.py -a $annotation -o $annotation_intron_list"
python gen_intronlist_from_annotation.py -a $annotation -o $annotation_intron_list
echo ""

echo "Compare spliced alignment to annotated introns:"
echo "python compare_intron_lists.py -O 20 -I $annotation_intron_list -C ${annotation_intron_list}.cov -a $alignment -M $max_intron_len -H $comparison_out -g $coverage_out"
python compare_intron_lists.py -O 20 -I $annotation_intron_list -C ${annotation_intron_list}.cov -a $alignment -M $max_intron_len -H $comparison_out -g $coverage_out

# Clean up
rm -r $tmpoutdir
