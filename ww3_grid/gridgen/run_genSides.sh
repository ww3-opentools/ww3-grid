#!/bin/bash

# run the grid cell face generation
# inputs are:
#  $1 working directory, the *Cels.dat and *.nml file(s) need to be in here
#  $2 the namelist to be used; e.g. smcGrid.nml, arcGrid.nml 

set -eu

WRKDIR=$1
NMLIST=$2

MYDIR=$PWD

cp genSides.exe $WRKDIR/work.genSides.exe
cp countijsdnew $WRKDIR/work.countijsdnew

cd $WRKDIR

CHKARC=`grep ARCTIC $NMLIST | cut -f 2 -d '.'`
if [ "$CHKARC" == "TRUE" ]; then
   LOGFILE='arcSides_output.txt'
else
   LOGFILE='smcSides_output.txt'
fi
echo "[INFO] genSides logfile is: $LOGFILE"

cp $NMLIST smcSides.nml
echo '[INFO] Launching genSides.exe to generate face arrays'
work.genSides.exe >$LOGFILE

echo '[INFO] Sorting the face arrays'
CHKARC=`grep ARCTIC $NMLIST | cut -f 2 -d '.'`
if [ "$CHKARC" == "TRUE" ]; then
   work.countijsdnew ww3AISide.d ww3AJSide.d >>$LOGFILE
else
   work.countijsdnew ww3GISide.d ww3GJSide.d >>$LOGFILE
fi

echo '[INFO] Tidying up'
rm work.genSides.exe
rm work.countijsdnew
rm smcSides.nml
rm *ISide.d *JSide.d

cd $MYDIR

exit 0
