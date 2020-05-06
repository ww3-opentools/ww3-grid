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
# run_smcProps.sh
#
# PURPOSE:
#  Runs smcProps.f90 executable, for the smc propagation test and move results data
#  into smcProps subdirectory
#  inputs are:
#   $1 working directory, the *Cels.dat *Side.dat and *.nml file(s) need to be here
#
# REVISION HISTORY:
#
# J.G. Li; Met Office; Aug-2007; Version:0.1
#  Initial functions development and testing at Met Office
#
# A. Saulter; Met Office; May-2020 Version:1.0
#  Code library prepared for initial release on github
#
#==================================================================================

set -eu

if [ $# -ne 1 ]
  then echo "$0: Usage is $0 PATHTOWORKINGDIRECTORY "
  exit 1
fi
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
mv CMesgs.txt smcProps
cd smcProps
echo '[INFO] Propagation test data file list' > cfiles.txt
ls Cn*.d >> cfiles.txt

cd $MYDIR

exit 0
