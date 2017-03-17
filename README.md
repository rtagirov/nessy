# PMOD version of the code

The state of the NESSY code as it was on 01.10.2016 on the PMOD servers.
This branch runs with the 13th version of the IFORT compiler.
It does not run with the latest, i.e. 17th, version of IFORT (ICL servers as of 01.10.2016).
The HMINUS part of the code, however, also runs with the 15th version of IFORT (MPS servers as of 01.10.2016).
No tests were performed for the FIOSS part of the code and the 15th version of IFORT.

Compilation:

The object files (*.o and *.mod) are collected in the directory ./obj (see ./src/make/Makefile.opt; OUT=../obj/)

The compilation is parallelized between 4 processors (see ./src/Makefile; OPTS="-j4")

The list of the source files for FIOSS and HMINUS is given in ./src/make/Makefile.sources, this file is edited manually

The IFORT compiler options are given in ./src/make/Makefile.opt

Before compiling the code make sure to update the list of dependencies by running bash ./src/make/write_deps.sh ./src/*.for > ./src/make/Makefile.deps

To compile HMINUS + FIOSS, run "make all" in the ./src directory
