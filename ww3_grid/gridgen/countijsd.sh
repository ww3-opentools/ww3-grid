#!/bin/bash

#==================================================================================
# BSD License
#
# Copyright (c)2010-2020, ww3-opentools developers, all rights reserved
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
# countijsd.sh
#
# PURPOSE:
#  Count ISD JSD face numbers of different sizes
#  Usage:  countisd ISD_File JSD_File
#  Output is saved as newISD_File, newJSD_File
#
# REVISION HISTORY:
#
# J.G. Li; Met Office; Feb-2010; Version:0.1
#  Initial code development and testing at Met Office
#
# A. Saulter; Met Office; May-2020; Version:1.0
#  Code prepared for initial release on github
#
#==================================================================================

echo ''
echo '  Cell sorting routine countijsd'

##  Use given input files.
#
if [ $# -ne 3 ]
  then echo "$0: Usage is $0 NLEVS U_iside v_jside "
  exit 1
fi
#
NLEVS=$1
echo "Sorting sides for model with ${NLEVS} levels"
#
USide=$2
if test ! -f "$USide"  # Check if USide file exists
  then
  echo "  $USide does not exist"
  exit 1
fi
#
VSide=$3
if test ! -f "$VSide"  # Check if VSide file exists
  then
  echo "  $VSide does not exist"
  exit 1
fi
echo "  U/V side arrays are from: $USide $VSide"

##  sort according to y-size then j and i count
sort -s -k 3,3n -k 2,2n -k 1,1n  $USide > temp1
sort -s -k 8,8n -k 2,2n -k 1,1n  $VSide > temp2 

##  cut out y-size field for counting
##  Use awk field number to extract y-size column.
awk '{print $3}' temp1 > temp3
awk '{print $8}' temp2 > temp8

##  count the different sizes
# U cells
CSIZE=1
LEVS=1
CSTR=`cat temp3 | wc -l `
while [ $LEVS -le $NLEVS ]
do
  NCELLS=`grep '^'${CSIZE} temp3 | wc -l `
  CSTR=$CSTR' '$NCELLS
  let LEVS=$LEVS+1
  let CSIZE=$CSIZE*2
done
echo 'U Cell Counts; Total, N levels:'$NLEVS
echo $CSTR
echo $CSTR > temp4

# V cells
CSIZE=1
LEVS=1
CSTR=`cat temp8 | wc -l `
while [ $LEVS -le $NLEVS ]
do
  NCELLS=`grep '^'${CSIZE} temp8 | wc -l `
  CSTR=$CSTR' '$NCELLS
  let LEVS=$LEVS+1
  let CSIZE=$CSIZE*2
done
echo 'V Cell Counts; Total, N levels:'$NLEVS
echo $CSTR
echo $CSTR > temp5

##  Merge saved counts with sorted face array files
echo '  Writing sorted I face arrays to '$USide'at'
cat temp4 temp1 > $USide'at'
echo '  Writing sorted J face arrays to '$VSide'at'
cat temp5 temp2 > $VSide'at'

##  Clear temporary files
rm  temp[1-5] temp8

##  All done.
echo '  countijsd completed'
exit 0

