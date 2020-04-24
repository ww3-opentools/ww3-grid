#!/bin/bash

# run the smc propagation test and move data into smcProps subdirectory
# inputs are:
#  $1 working directory, the *Cels.dat *Side.dat and *.nml file(s) need to be in here

set -eu

WRKDIR=$1

MYDIR=$PWD

cp smcProps.exe $WRKDIR/work.smcProps.exe

echo "[INFO] Running propagation tests in $WRKDIR"
cd $WRKDIR
echo '[INFO] Launching smcProps.exe for propagation tests'
echo "[INFO] smcProps logfile is: smcProps_output.txt"
work.smcProps.exe >smcProps_output.txt

echo '[INFO] Tidying up'
echo "[INFO] Moving smcProps files to ${WRKDIR}/smcProps"
if [ ! -d smcProps ]; then
  mkdir smcProps
fi
rm work.smcProps.exe
mv smcProps_output.txt smcProps
mv Cn*.d smcProps
cd smcProps
echo '[INFO] Propagation test data file list' > cfiles.txt
ls Cn*.d >> cfiles.txt

cd $MYDIR

exit 0
