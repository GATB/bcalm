#!/bin/bash
    
rm -f reference.fasta
wget https://raw.githubusercontent.com/GATB/MindTheGap/f5cb0fec816686c7393772787d736565c4f056a4/test/full_test/reference.fasta >& /dev/null
../build/bcalm -in reference.fasta -abundance-min 1 > /dev/null >& /dev/null
rm -rf reference.unitigs.fa.glue*
rm -f compare_fasta.py
wget https://raw.githubusercontent.com/GATB/minia/master/test/compare_fasta.py >& /dev/null
python compare_fasta.py reference.fasta reference.unitigs.fa
res=$?
rm -f reference.fasta reference.h5 reference.unitigs.fa compare_fasta.py 
if [ "$res" = "0" ]
then
    echo "test OK"
    exit 0
else
    echo "test KO"
    exit 1
fi
