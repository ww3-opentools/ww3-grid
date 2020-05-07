#!/bin/bash

#==================================================================================
# BSD License
#
# Copyright (c)2007-2020, ww3-opentools developers, all rights reserved
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice, this
#  list of conditions and the following disclaimer in the documentation and/or
#  other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#  contributors may be used to endorse or promote products derived from this
#  software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
#
#==================================================================================
# run_genSides.sh
#
# PURPOSE:
#  Runs genSides.f90 executable and countijsdnew scripts to generate SMC grid cell
#  face arrays for use with WAVEWATCH III
#
# REVISION HISTORY:
#
# J.G. Li; Met Office; Aug-2007
#  Initial functions development and testing at Met Office
#
# A. Saulter; Met Office; May-2020
#  Code library prepared for initial release on github
#
#==================================================================================

set -eu

# inputs are:
#  $1 working directory, the *Cels.dat and *.nml file(s) need to be in here
#  $2 the namelist to be used; e.g. smcGrid.nml, arcGrid.nml 
WRKDIR=$1
NMLIST=$2

MYDIR=$PWD

cp genSides $WRKDIR/work.genSides
cp countijsd.sh $WRKDIR/work.countijsd.sh

cd $WRKDIR

CHKARC=`grep ARCTIC $NMLIST | cut -f 2 -d '.'`
if [ "$CHKARC" == "TRUE" ]; then
   LOGFILE='arcSides_output.txt'
else
   LOGFILE='smcSides_output.txt'
fi
echo "[INFO] genSides logfile is: $LOGFILE"

cp $NMLIST smcSides.nml
echo '[INFO] Launching genSides to generate face arrays'
work.genSides >$LOGFILE

echo '[INFO] Sorting the face arrays using countijsd'
NLEVS=`grep NLEVS $NMLIST | cut -f 2 -d '='`
CHKARC=`grep ARCTIC $NMLIST | cut -f 2 -d '.'`
if [ "$CHKARC" == "TRUE" ]; then
   work.countijsd.sh $NLEVS ww3AISide.d ww3AJSide.d >>$LOGFILE
else
   work.countijsd.sh $NLEVS ww3GISide.d ww3GJSide.d >>$LOGFILE
fi

echo '[INFO] Tidying up'
rm work.genSides
rm work.countijsd.sh
rm smcSides.nml
rm *ISide.d *JSide.d

cd $MYDIR

exit 0
