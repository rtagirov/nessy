# ICL version

The state of the code on 25.12.2016, when it started running on my machine at ICL
(IFORT version 17.0.0 20160721).
Below 70 nm the resulting spectrum produced by this version yields a discrepancy of about 8%
with the pmod branch spectrum.

Compilation:

The object files (*.o and *.mod) are collected in the directory ./obj (see ./src/make/Makefile.opt; OUT=../obj/)

The compilation is parallelized between 16 processors (see ./src/Makefile; OPTS="-j16")

The list of the source files for FIOSS and HMINUS is given in ./src/make/Makefile.sources, this file is edited manually

The IFORT compiler options are given in ./src/make/Makefile.opt

Before compiling the code make sure to update the list of dependencies by running bash ./src/make/write_deps.sh ./src/*.for > ./src/make/Makefile.deps

To compile HMINUS, FIOSS or HMINUS + FIOSS, run, respectively, the 3 following commands in the ./src directory:

1. make hminus

2. make fioss

3. make all

Note: the 32-bit architecture option has been removed from ./src/Makefile and ./src/make/Makefile.opt
