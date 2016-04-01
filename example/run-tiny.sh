../build/bcalm -in tiny_read.fa -k 12 -abundance 1 $1 $2
../build/bglue -in unitigs.h5 -k 12 $1 $2
