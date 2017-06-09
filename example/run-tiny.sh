bcalm=`find ../*/bcalm | head -n 1` # covers bin/ and build/ cases
$bcalm -in tiny_read.fa -kmer-size 13 -abundance-min 1 $1 $2
