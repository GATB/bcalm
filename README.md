# Welcome

This project has been automically created with a gatb-core script.

It can be used as a starter for developping tools based on gatb-core.

The architecture of the project is:

    * a CMakeLists.txt file used for building the project
    * a 'src' file holding a default main function
    * a 'thirdparty' directory holding the gatb-core resources


# Project build

For building your project, you should do the following
    
    mkdir build;  cd build;  cmake ..;  make
    
Then, you should get a binary holding the name of the project.

Note: the first compilation should take some time since the GATB-CORE library is generated.



# Documentation

If doxygen is installed, you can generate the gatb-core documentation with 'make doc'

The documention is available in HTML, the entry point being 'ext/gatb-core/doc/html/index.html'



# Examples

The project is created with a default 'main' function that dumps some information about the library.

You can find many snippets showing how to use the library. 
These snippets are located in 'thirdparty/gatb-core/examples' and are split in several fields.

You can copy the content of one of the snippet file into the 'src/main.cpp' file and relaunch the build.
For instance:
    cp thirdparty/gatb-core/examples/debruijn/debruijn4.cpp src/main.cpp

WARNING... Some examples use on purpose lambda expressions, so you will need a compiler supporting this feature for this examples.


# Binaries from gatb-core

After the project build, some gatb-core binaries are available here 'ext/gatb-core/bin'

As gatb-core uses HDF5, you will have here some H5xxx tools from the HDF5 distribution.

You will also find two gatb-core binaries:

    * dbgh5:    builds a DeBruijn graph from a set of reads and save it as a HDF5 file
    * dbginfo:  dumps information about a DeBruinj graph build by dbgh5

