#!/bin/bash

set -e

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % SAFT examples: get_data.sh %
echo % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo 
echo This script gets reads and the genome sequence 
echo for two C. elegans examples.
echo 

export DATA_DIR=data

if [ -d $DATA_DIR ]
then
	echo Data directory ./$DATA_DIR already exists
	echo \(remove it to run this script again\)
	exit -1
fi

echo Downloading rQuant example data from FTP server ...
wget -c ftp://ftp.tuebingen.mpg.de/fml/group_raetsch/software/SAFT/SAFT_examples.tar.gz
echo uncompressing ...
tar -xzf SAFT_examples.tar.gz 
mv SAFT_examples $DATA_DIR
echo
echo -n Done.
echo
