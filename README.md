# LTE Marking

Note: currently this branch is almost a copy of the acc_damp branch, its development is underway

The following changes with respect to the master branch have been implemented:

1. The acceleration damping, which is necessary to run the code for certain types of atmospheres (see comblock.for, etl.for, como.for, phys.for; search for
   "damp_acc", "acc_damp", "opt_dep"); to switch off the acceleration damping use the variable "damp_acc" in comblock.for

2. The Rayleigh scattering on neutral hydrogen (see comblock.for, opac.for, phys.for; search for "rayleigh" or "sigma_rayleigh"); to switch off the rayleigh scattering
   use the variable "rayleigh" in comblock.for

3. The FAL_VD and TABLE files are abandoned. The FAL and TABLE options in the CARDS file are obsolete now. Any model of the atmosphere should be named ATM_MOD.
   The code automatically understands what format of model it is (FAL, Kurucz or MURAM) by the number of columns the atmosphere model file has
   (7 - Kurucz, 5 - FAL, 4 - MURAM). The code reads the file and prints out the atmospheric stratification file ATM_STR containing various quantities of interest
   pertaining to the atmopheric parameters stratification (gradients, height grid step, etc.).
   See geomesh.for, fileoper.for; search for "print_strat", "read_atm_file_col", "read_fal_fmt_mod", "read_kur_fmt_mod", "read_mur_fmt_mod".

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

Note: if "rayleigh" and "damp_acc" variables are set to .false. the functioning of this version is the same as that of the master branch
