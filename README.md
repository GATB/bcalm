# BCALM 2 

BCALM 2 is a bioinformatics software for constructing the compacted de Bruijn graph of sequencing data.

This repository is the new, parallel version of the BCALM software.
It is using a new algorithm, and is implemented using the GATB library. 
The original, single-threaded code of BCALM (version 1) is still available at: https://github.com/Malfoy/bcalm

[![Build Status](https://travis-ci.org/GATB/bcalm.svg?branch=master)](https://travis-ci.org/GATB/bcalm)

# Usage

Read the instructions below to compile, then:

    ./bcalm -in [reads.fa] -kmer-size [kmer_size] -abundance-min [abundance_threshold]
  
e.g.

    ./bcalm -in reads.fastq -kmer-size 21 -abundance-min 2

Importants parameters are:

    -kmer-size [int]
    
The k-mer size, i.e. length of the nodes of the de Bruijn graph.

    -abundance-min [int]

Sets a threshold X below which k-mers that are seen (strictly) less than X times in the dataset are filtered out; i.e. sequencing errors, typically.

# Pre-requisites:

GCC >= 4.8 or a very recent C++11 capable compiler

# Installation

To retrieve bcalm and its submodule (gatb-core), type

    git clone --recursive https://github.com/GATB/bcalm
    
# Input formats

File input format can be fasta, fastq, either gzipped or not.

To pass several files as input, separate file names by a comma (","), for example:

    ./bcalm -in A1.fa,A2.fa,A3.fa [..]

Alternatively, input can be a list of files (one file per line):

    ls -1 *.fastq > list_reads
    ./bcalm -in list_reads [..]
   
# Output

BCALM 2 outputs the set of unitigs of the de Bruijn graph.
A unitig is the sequence of a non-branching path. Unitigs that are connected by an edge in the graph overlap by exactly (k-1) nucleotides. 
FASTA output header: 

    ><id> LN:i:<length> KC:i:<abundance> KM:f:<abundance> L:<+/->:<other id>:<+/-> [..]

Where:

* `LN` field is the length of the unitig
    
* `KC` and `KM` fields are for total abundance and mean abundance of kmers inside the unitig, respectively.

* Edges between unitigs are reported as `L:x:y:z` entries in the FASTA header (1 entry per edge). A classic forward-forward outcoming edge is labeled `L:+:[next node]:+`. A forward-reverse, `L:+:[next node]:-`. Incoming edges are labeled as follows: `L:-:[previous node]:+`, meaning that if you reverse-complemented the current node, then there would be an edge from the last k-mer of current node to the first k-mer of the forward strand of [previous node]. (Think of incoming edges as outcoming edges of the reverse-complement of the node.)
 
# Reverse-complements

BCALM 2 converts all k-mers into their canonical representation with respect to reverse-complements.
In other words, a k-mer and its reverse complement are considered to be the same object, appearing only once in the output, either in forward or reverse orientation.

# Project build

For building your project, you should do the following
    
    mkdir build;  cd build;  cmake ..;  make -j 8
    
Then, you should get a binary in

    tools/bcalm

Note: the first compilation should take some time since the GATB-CORE library is generated.

# Larger k values

BCALM 2 supports arbitrary large k-mer lengths. You need to recompile it from sources. For k up to, say, 320, type this in the build folder:

    rm -Rf CMake* && cmake -DKSIZE_LIST="32 64 96 128 160 192 224 256 320" .. && make -j 8

The list of kmers should only contain multiples of 32. Intermediate values create optimized code for smaller $k$'s. You could specify just `KSIZE_LIST="320"` but then using smaller k values would be as slow as large ones.

Acknowledgements
========
If using BCALM 2, please cite:
Rayan Chikhi, Antoine Limasset and Paul Medvedev, Compacting de Bruijn graphs from sequencing data quickly and in low memory, Proceedings of ISMB 2016, Bioinformatics, 32 (12): i201-i208. 

This project has been supported in part by NSF awards DBI-1356529, CCF-1439057, IIS-1453527, and IIS-1421908.
