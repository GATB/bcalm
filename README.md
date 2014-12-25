# Welcome

This is a port of BCALM using the GATB library. Under development.

# Installation

To retrieve bcalm and its submodule (gatb-core), type

    git clone --recursive https://github.com/GATB/bcalm

# Project build

For building your project, you should do the following
    
    mkdir build;  cd build;  cmake ..;  make
    
Then, you should get a binary in

    tools/bcalm

Note: the first compilation should take some time since the GATB-CORE library is generated.


# Documentation

If doxygen is installed, you can generate the gatb-core documentation with 'make doc'

The documention is available in HTML, the entry point being 'ext/gatb-core/doc/html/index.html'


