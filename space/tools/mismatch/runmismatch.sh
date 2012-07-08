#!/bin/bash
cd /fml/ag-raetsch/nobackup/projects/rgasp/reads
for i in `find ./ -name read_sample.result` ; do python ~/svn/projects/rgasp/read_processing/mismatch/mismatch.py $i `dirname $i`/mismatch.pickle ; done
for i in `find ./ -name mismatch.pickle` ; do python ~/svn/projects/rgasp/read_processing/mismatch/plot_mismatch.py $i `dirname $i`/mismatch.pdf ; done
