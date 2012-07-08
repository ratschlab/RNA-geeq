#!/bin/sh

set -e

# IN-FILES
if [ -z "$1" ]
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<statistics_outfile\>
    exit 2
else
    alignment=$1
fi

if [ -z "$2" ]
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<statistics_outfile\>
    exit 2
else
    genome_fasta=$2
fi

# OUT-FILES
if [ -z "$3" ]
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<statistics_outfile\>
    exit 2
else
    statistics_out=$3
fi

# INTERMEDIATE-FILES
tmpoutdir=`mktemp -d` || exit 1
annotation_intron_list=$tmpoutdir/annotation.introns

# Check for genome_config
if [ ! -d "$tmpoutdir/genome" ];
then
    echo "creating genome information object from $genome_fasta in $tmpoutdir ..."
        python ./tools/make_gio.py $genome_fasta $tmpoutdir
    echo "... done"
    echo ""
fi

genome_config="$tmpoutdir/genome.config"

echo "Generate statistics:"
echo "python gen_alignment_statistics.py -a $alignment -g $genome_config -H $statistics_out"
python gen_alignment_statistics.py -a $alignment -g $genome_config -H $statistics_out
echo ""

# Clean up
rm -r $tmpoutdir
