#!/bin/sh

set -e

# IN-FILES
if [ -z "$1" ]
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<output_file\>
    exit 2
else
    alignment=$1
fi

if [ -z "$2" ]
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<output_file\>
    exit 2
else
    genome_fasta=$2
fi

# OUT-FILES
if [ -z "$3" ]
then
    echo usage $0 \<alignment\> \<genome_fasta\> \<output_file\>
    exit 2
else
    outfile=$3
fi

# INTERMEDIATE-FILES
tmpoutdir=`mktemp -d` || exit 1

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
echo "python evaluateAlignmentError.py -a $alignment -g $genome_config -H $outfile"
python evaluateAlignmentError.py -a $alignment -g $genome_config -H $outfile
echo ""

# Clean up
rm -r $tmpoutdir
