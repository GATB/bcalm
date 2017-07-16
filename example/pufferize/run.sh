../../build/bcalm -in refs.fa -abundance-min 1 -kmer-size 9 -minimizer-size 5
python ../../scripts/pufferize.py refs.fa refs.unitigs.fa 9
