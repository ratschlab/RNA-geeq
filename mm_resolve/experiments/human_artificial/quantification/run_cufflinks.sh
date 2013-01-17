#!/bin/bash

set -e

usage() {
    echo "usage: $0 <genes> <readnum> <experiment> <stage> [noise] [<chr1,[chr2,...]>]"
    exit 1
}

[[ -z "$1" ]] && usage
genes="$1"
shift

[[ -z "$1" ]] && usage
readnum="$1"
shift

[[ -z "$1" ]] && usage
which_set="$1"
shift

[[ -z "$1" ]] && usage
stage="$1"
shift

noise="$1"
noise_out=""
if [ ! -z "$noise" ]
then
    noise_out=".noise$noise"
    noise="noise${noise}."
fi

chrms="$2"
if [ ! -z "$chrms" ]
then
    chrms=`echo $chrms | tr ',' '_'`_
fi

work_dir=${genes}_genes_${readnum}_reads
out_dir=${work_dir}/cufflinks/${which_set}${noise_out}.stage${stage}

gtf=`pwd`/annotation/hg19_${chrms}subsample_${genes}_genes.nochr.gtf

if [ "${which_set}" != "unfiltered" ]
then
    which_set="${which_set}."
else
    which_set=""
fi
alignment=`pwd`/${genes}_genes_${readnum}_reads/hg19_${chrms}subsample_${genes}_genes.gtf.${noise}fastq.gz.mapped.${stage}.${which_set}sorted.bam

export LD_LIBRARY_PATH=/cbio/grlab/share/software/lib/:$LD_LIBRARY_PATH

mkdir -p ${out_dir}
if [ ! -f ${out_dir}/cufflinks.log ]
then
    /cbio/grlab/share/software/cufflinks/cufflinks-1.3.0.Linux_x86_64/cufflinks -o ${out_dir} -G $gtf ${alignment} 2> ${out_dir}/cufflinks.log 
    python cuff_extract_counts.py ${out_dir}/transcripts.gtf > ${out_dir}/transcripts.counts
else
    echo "${out_dir}/cufflinks.log exists - delete to re-run"
fi

