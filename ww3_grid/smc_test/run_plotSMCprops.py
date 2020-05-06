#==================================================================================
# BSD License
#
# Copyright (c)2018-2020, ww3-opentools developers, all rights reserved
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
# run_plotSMCprops.py
#
# PURPOSE:
#  Runs plotSMCgrid.py functions, based on namelist file inputs, for plotting 
#  results of SMC grid propagation tests
#
# REVISION HISTORY:
#
# J.G. Li; Met Office; Dec-2018; Version:0.1
#  Initial functions development and testing at Met Office
#
# A. Saulter; Met Office; May-2020; Version:1.0
#  Code library prepared for initial release on github
#
#==================================================================================

import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from datetime import datetime, timedelta
import sys

from plotSMCgrid import readtext
from plotSMCgrid import readcell   
from plotSMCgrid import steromap
from plotSMCgrid import rgbcolor
from plotSMCgrid import swhglobl
from plotSMCgrid import swhlocal

def readGridnml(Wrkdir='.',fname='smcGrid.nml'):
    """Read grid information from smcGrid namelist file"""

    gridmeta={}
    nmlfile = Wrkdir+'/'+fname
    runme=True
    if runme:
    #try:
        with open(nmlfile,'r') as inp:
            rdfile = inp.readlines()
            inp.close()

        gridmeta['nlevs'] = np.int(rdfile[1].split('=')[1])
        gridmeta['dx'] = np.float(rdfile[4].split('=')[1])
        gridmeta['dy'] = np.float(rdfile[5].split('=')[1]) 
        gridmeta['x0'] = np.float(rdfile[6].split('=')[1])
        gridmeta['y0'] = np.float(rdfile[7].split('=')[1]) 
        gridmeta['nx'] = np.int(rdfile[3].split('=')[1])
        gridmeta['ny'] = np.int(rdfile[2].split('=')[1]) 
        gridmeta['grdfile'] = rdfile[9].split('=')[1].split("'")[1] 
        # use RUNARCTIC to tell me if I need to load an arctic grid file
        if rdfile[18].split('=')[1].split('.')[1].lower() == 'true':
            gridmeta['arcfile'] = rdfile[19].split('=')[1].split("'")[1] 
        else:
            gridmeta['arcfile'] = None
        print('[INFO] Read grid data from %s' %nmlfile)
    #except:
    #    print('[ERROR] Unable to find file %s' %nmlfile)

    return gridmeta


def pltProps(Wrkdir, ModlName, Cel_file='ww3Cels.dat', 
              Arc_file=None, DatGMC=None, Flsdir=None, figtype='.png'):

    """Plots results of SMC propagation tests"""
    # Note - postscript is presently hardwired in swhglobl and swhlocal

    if DatGMC is None:
        DatGMC=Wrkdir
    if Flsdir is None:
        Flsdir=Wrkdir+'/smcProps'
    Cel_file = Wrkdir + '/' + Cel_file
    if Arc_file is not None:
        Arc_file = DatGMC+'/'+Arc_file

    if Arc_file is not None:
        headrs, cel = readcell( [Cel_file, Arc_file] ) 
        na = int( headrs[1].split()[0] )
        nb = int( headrs[1].split()[1] )
    else:
        headrs, cel = readcell( [Cel_file] ) 
        na = 0
        nb = 0
    ng = int( headrs[0].split()[0] )
    nc = ng + na
    print('[INFO] Merged total cel number = %d' % nc )

    # Maximum j row number in Global part
    jmxglb = cel[ng-1,1] 
    print ('[INFO] Maximum j row in global grid = %d' % jmxglb )

    # Use own color map and defined depth colors 
    colrfile = 'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

    ##  Possible selection of your plot types. 
    gorloc={0:'Global',1:'EuroArc',2:'Pacific',3:'Atlantic'}

    # Prompt selection choices and ask for one input
    print(" \n ", gorloc)
    instr = input("\n *** Please enter your selected number here > ")
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print("Invalid selection, program terminated.")
        exit()

    print(" Draw SWH plots "+pltype)

    # Choose global or local verts from different files.
    vrfile = DatGMC+'/'+ModlName+'_Vrts'+pltype[0:4]+'.npz'
    vrtcls = np.load( vrfile )

    if( pltype == 'Global' ):
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']
        svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
        config = vrtcls['cnfg']
        print('[INFO] n/svrts/cels config read ')
    else:
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel'] 
        config = vrtcls['cnfg']
        print('[INFO] nvrts, ncels and config read ')

    # EuroArc and Pacific paper orientations
    if( pltype == 'Pacific' ):
        papror='landscape'
    else:
        papror='portrait'

    # Define spectral direction
    ndir=36
    theta=np.arange(ndir)*np.pi*2.0/ndir

    # Add a spectral array plots for the Northern stripe 
    x0= 3.0
    y0= 5.0
    t0=np.pi*0.25
    cs=np.cos(theta + t0)
    xn=theta*0.0+x0
    yn=theta*0.0+y0
    for i in range(ndir):
        if( cs[i] > 0.0 ): 
            spc=1.2*cs[i]*cs[i]
            xn[i]=x0+spc*np.cos(theta[i])
            yn[i]=y0+spc*np.sin(theta[i])

    # Add another spectral array plots for the Southern stripe 
    x1=-2.0
    y1=-8.0
    cs=np.cos(theta - t0)
    xs=theta*0.0+x1
    ys=theta*0.0+y1
    for i in range(ndir):
        if( cs[i] > 0.0 ): 
            spc=1.2*cs[i]*cs[i]
            xs[i]=x1+spc*np.cos(theta[i])
            ys[i]=y1+spc*np.sin(theta[i])

    #  Polar disk spectral array uses the Southern one but new location
    xp= 3.0
    yp= 8.0

    #  Atlantic disk spectral array uses the Southern one but new location
    xt=-1.0
    yt=-2.0

    # Use ijk to count how many times to draw.
    ijk=0

    # Specify number of steps per hour
    nhr=24

    # Read in cell concentration data files from a list file
    hdr, cnfiles = readtext(Flsdir+'/cfiles.txt')
    cfiles = cnfiles.astype(np.str).reshape(len(cnfiles))

    # loop over available files 
    for nn in range(0,len(cnfiles),1):
        dfile=Flsdir+'/'+cfiles[nn] 

        hdlist, swh2d = readtext(dfile)
        mt = int(hdlist[0])
        mc = int(hdlist[1])
        swhs = swh2d.flatten()[0:mc]

        # Skip Arctic polar cell if nc = nga
        if( mc != nc ):
            print( '[INFO] Unmatching mc/nc = %d %d' % (mc, nc) ) 
            exit()
        else:
            print('[INFO] Plotting cell number mc = %d' % mc )

        # Convert time step for output file
        ntsp='NTS = %5d' % (mt)
        thrs='T = %d hours' % (np.int(cfiles[nn][2:7])-10000) 

        # Call function to draw the swh plot.
        figfl = Flsdir + '/Hs' + cfiles[nn][2:7] + figtype
        if( pltype == 'Global' ):
            swhglobl(swhs, nvrts, ncels, svrts, scels, colrs, config,
                     mdlname= ModlName, datx=thrs, psfile=figfl)
        else:
            swhlocal(swhs,nvrts, ncels, colrs, config,
                     mdlname= ModlName, datx=thrs,psfile=figfl,
                     paprorn=papror )

    # Increase ijk for next plot
        ijk += 1
        print("[INFO] Finish plot No.", ijk," at ", datetime.now())

    # End of date loop
    return


#---
# main program

ModlName = sys.argv[1]
Wrkdir   = sys.argv[2]

# read the grid metadata from the grid namelist file
gridmeta = readGridnml(Wrkdir,fname='smcGrid.nml')

# plot the propagation test data
pltProps(Wrkdir, ModlName, 
         Cel_file=gridmeta['grdfile'], 
         Arc_file=gridmeta['arcfile'], figtype='.ps')
