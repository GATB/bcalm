# Welcome

BCALM is a bioinformatics software for constructing the de Bruijn graph of sequencing data.

More precisely, this repository is the next version of the BCALM software (original code: https://github.com/Malfoy/bcalm) using the GATB library. It is under development.

# Installation

To retrieve bcalm and its submodule (gatb-core), type

    git clone --recursive https://github.com/GATB/bcalm

# Project build

For building your project, you should do the following
    
    mkdir build;  cd build;  cmake ..;  make -j 8
    
Then, you should get a binary in

    tools/bcalm

Note: the first compilation should take some time since the GATB-CORE library is generated.
