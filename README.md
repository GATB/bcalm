# Welcome

BCALM is a bioinformatics software for constructing the de Bruijn graph of sequencing data.

More precisely, this repository is the next version of the BCALM software (original code: https://github.com/Malfoy/bcalm) using the GATB library. 

# Under development

Need to check correctness. Output might be missing a few unitigs due to edge cases.

# Usage

Read the instructions below to compile, then:

    ./bcalm -in [reads.fa] -k [kmer_size] -abundance [abundance_threshold]
    ./bglue -in unitigs.h5 -k [kmer_size]
  
e.g.

    ./bcalm -in reads.fa -k 21 -abundance 2
    ./bglue -in unitigs.h5 -k 21
 

# Installation

To retrieve bcalm and its submodule (gatb-core), type

    git clone --recursive https://github.com/GATB/bcalm

# Project build

For building your project, you should do the following
    
    mkdir build;  cd build;  cmake ..;  make -j 8
    
Then, you should get a binary in

    tools/bcalm

Note: the first compilation should take some time since the GATB-CORE library is generated.

