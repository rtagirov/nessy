#!/bin/bash
#
# Generate the dependencies for the given files using modules
# for every file A, search for use mod_*, strip of the use mod_ part and write out A: $((addprefix) $(OUT), "dependencies" )
#
[[ "$@" == "" ]] && echo "Usage: $0 file1 ... filen > Makefile.deps

with fileN the fortran files to check for modules.

E.g. $0 *.for

The programm looks for the lines use MOD_..., strips off the MOD_ part and writes out a
fileN: \$(addprefix \$(OUT), MODULES ).
" || MOD_FILE=`mktemp -t write_deps.XXXXX`

## Pretty print...
MAX_LENGTH=0
for  A in $@; do [[ $(echo $A | wc -L)  -gt MAX_LENGTH ]] && MAX_LENGTH=$(echo $A | wc -L); done
MAX_LENGTH=$[MAX_LENGTH+1]
SPACE="                                                                                                      "
SPACE=${SPACE:1:$MAX_LENGTH}




for A in $@; do
  # grep the module name, remove module and any comments, etc
  echo $(echo $A | sed 's/\.for$/.o/'):$(grep -i '^ *module ' $A | sed -e 's/^ *module //i;s/[,!].*//i;s/^ *//;s/ //i;s/[^A-Za-z0-9_]//i'|tr [:upper:] [:lower:]) >> $MOD_FILE
done

#go through all the files
for A in $@; do
  # go through every USE statement
  MODS=""
  for B in $(grep -i '^ *use ' $A | sed -e 's/^ *use //i;s/[,!].*//i;s/[^A-Za-z0-9_]//i'); do
    # search the MOD_FILE for the modules in the use statement
    MODS="$MODS $(grep -i ":$B"  $MOD_FILE  | sed -e 's/:.*//;s/ //;s/\n//')"
  done
#echo $A::$MODS
echo $MODS | grep '^ *$' &> /dev/null || echo "$A: ${SPACE:${#A}}  \$(addprefix \$(OUT)," $MODS ")"

done
rm $MOD_FILE
