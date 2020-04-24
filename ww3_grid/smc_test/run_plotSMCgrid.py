"""
##
##  It reads cell files and uses own projection to draw 
##  the SMC grid. Projected polygons are collected into 
##  vert variables for grid and subsequent swh plots. 
##
#;  First created: For 4 views     J G Li   26 Nov 2008
#;  Modified for Color grids       J G Li   26 Aug 2009
#;  Adapted for 25km SMC grids     J G Li    3 Feb 2010
#;  Updated G25SMC-12c grids.      J G Li   25 May 2011
#;  Sterographic projection nR.    J G Li   30 Jun 2011
#;  Adapted for G50SMC grid.       J G Li   19 Aug 2011
#;  Extended to fill the Arctic.   J G Li    5 Oct 2011
#;  Rectify polar cell position.   J G Li   25 Oct 2011
#;  Simplify with readcell and steromap.  JGLi12Nov2014
##
##  Converted into a Python function.     JGLi05Dec2018
##  Save ELat/Lon and sx/yc in file.      JGLi11Dec2018
##  Add color map and draw color bar.     JGLi14Dec2018
##  Adapted for SMC36125 grid plot.       JGLi03Jan2019
##  Import resource and set stacksize.    JGLi07Jan2019
##  Use polycollections for two plots.    JGLi30Jan2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Adapted for Andy's SMC6125  grid.     JGLi12Mar2019
##  Adapted for generic SMC grids         ASaulter30Apr2020
##
"""

##  Import relevant modules and functions

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from datetime import datetime
import sys

from plotSMCgrid import readcell   
from plotSMCgrid import steromap
from plotSMCgrid import rgbcolor
from plotSMCgrid import smcglobl
from plotSMCgrid import smclocal


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


def pltSMCgrid(ModlName, nlevs, dx, dy, x0, y0, nx, ny,
                Wrkdir, Cel_file='ww3Cels.dat', Arc_file=None, DatGMC=None,
                figtype='.ps'):

    """Plots SMC grid data as postscript format file"""
    # Note - postscript is presently hardwired in smcglobl and smclocal

    print("[INFO] Program started at %s " % datetime.now().strftime('%F %H:%M:%S'))

    #  Read global and Arctic part cells. 
    if DatGMC is None:
        DatGMC=Wrkdir
    Cel_file = DatGMC+'/'+Cel_file
    if Arc_file is not None:
        Arc_file = DatGMC+'/'+Cel_file

    if Arc_file is not None:
        headrs, cel = readcell( [Cel_file, Arc_file] ) 
        na = int( headrs[1].split()[0] )
        nb = int( headrs[1].split()[1] )
    else:
        headrs, cel = readcell([Cel_file]) 
        na = 0
        nb = 0
    ng = int( headrs[0].split()[0] )
    nc = ng + na
    print('[INFO] Merged total cel number = %d' %nc)

    #  Size-1 cell increments
    dxlon = dx / 2.0**(nlevs-1)
    dylat = dy / 2.0**(nlevs-1)

    #  Origin and cells numbers for size-1 cell i j indices.
    xlon0 = x0 - dx/2.0
    ylat0 = y0 - dy/2.0
    nx = (nx+1) * 2**(nlevs-1) + 1
    ny = (ny+1) * 2**(nlevs-1) + 1

    #  Half the latitude grid to set j=0 on the Equator.
    neqt=0
    xlon=(np.arange(nx))*dxlon + xlon0
    ylat=(np.arange(ny))*dylat+ylat0
    print('[INFO] Full and Half lat grid = %d %d' %(ny-1, neqt))

    #  Maximum j row number in Global part
    imx = cel[:,0].max()
    jmx = cel[:,1].max()
    dimx= cel[:,2].max()
    djmx= cel[:,3].max()
    khmx= cel[:,4].max()
    print('[INFO] Maximum i j di dj k = %d %d %d %d %d' %(imx, jmx, dimx, djmx, khmx))
    khmn= cel[:,3].min()
    print('[INFO] Minimum depth in cel = %.1f' %khmn)

    #  Extra array in config for Arctic part
    ngabjm = [ng, na, nb, jmx]
    if Arc_file is None:
        Arctic = False

    #  Use own color map and defined depth colors 
    colrfile = 'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

    #  Maximum mapping radius.
    radius=10.0

    #  Possible selection of your plot types. 
    gorloc={0:'Global',1:'EuroArc',2:'Pacific',3:'Atlantic'}

    #  Prompt selection choices and ask for one input
    print(gorloc)
    instr = input(' *** Please enter your selected number here > ')
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print("Invalid selection, program terminated.")
        exit()

    print("[INFO] Draw SMC grid "+pltype)

    if( pltype == 'Global'):
        #  Whole global projection angle from N Pole to be 90.0
        pangle=90.0
        plon= 0.0 
        plat= 23.5 
        clrbxy=[ -9.6,-12.6, 19.0,  1.0]
        sztpxy=[ 16.0, 10.0,-10.2,-10.3]
        rngsxy=[-10.0, 10.0,-13.0, 10.0]

    if( pltype == 'EuroArc'):
        #  Euro-Arctic regional plot
        pangle=27.5 
        plon=  4.0
        plat= 69.0
        clrbxy=[  7.5, -4.5,   0.8,  7.0]
        sztpxy=[ 11.0, 15.0,   4.0, -5.0]
        rngsxy=[-10.0, 10.0, -14.0, 14.0]
        papror='portrait'

    if( pltype == 'Pacific'):
        #  West Pacific regional plot
        plon= 138.0
        plat=  18.0
        pangle=33.0 
        clrbxy=[-13.6,  7.5,   7.5,  0.8]
        sztpxy=[ 15.0, 10.0, -10.0,  7.0]
        rngsxy=[-15.0, 15.0, -10.0, 10.0]
        papror='landscape'

    if( pltype == 'Atlantic'):
        #  Whole global projection angle from N Pole to be 90.0
        pangle=64.0
        plon= 320.0 
        plat= 23.6 
        clrbxy=[ -9.6, -9.9, 19.0,  0.8]
        sztpxy=[ 10.0, 12.0, -6.4, -4.6]
        rngsxy=[-10.0, 10.0,-10.0, 11.0]
        papror='portrait'

    #print " Start loop over cells ... "
    print( "[INFO] Start loop over cells at %s " % datetime.now().strftime('%F %H:%M:%S') )

    #  Initial verts and ncels variable for polycollections.
    nvrts = []
    ncels = []
    svrts = []
    scels = []
    #  Exclude last cell, the North Polar Cell, and Arctic boundary cells.
    for i in range(nc): 
        if( (i < ng-nb) or (i >= ng+nb) ):
            xc=[cel[i,0],cel[i,0]+cel[i,2],cel[i,0]+cel[i,2],cel[i,0]]
            yc=[cel[i,1],cel[i,1],cel[i,3]+cel[i,1],cel[i,3]+cel[i,1]]
            slat=ylat[yc]
            slon=xlon[xc]

            #  Convert slat slon to elat elon with given new pole
            elat,elon,sxc,syc = steromap(slat,slon,plat,plon,Pangl=pangle,Onecl=True)

            if( (elat[0] >= 0.0) and (rngsxy[0] < sxc[0] < rngsxy[1])
                                 and (rngsxy[2] < syc[0] < rngsxy[3]) ):
                nvrts.append(list(zip(sxc,syc)))
                ncels.append(i)

            if( (elat[0] <  0.0) and (pltype == 'Global') ):
                svrts.append(list(zip(sxc,syc)))
                scels.append(i)

    #  End of cell i loop excluding polar cell
    print( "[INFO] Processing polar cell at %s " % datetime.now().strftime('%F %H:%M:%S') )
    #;  Polar cell as a octagon, note polar cell size set as size-64 for flux calculation
    #;  As it mergers 8 size-64 cells, each side of the octagon should be size-64 
    #;  10 apexes are calculated to conform with other cells.  To plot it, drop the last apex.
    #  Use a square box for polar cell to fit the new array size.  JGLi30Jan2019
    #i=nc-1
    #xc=cel[i,0]+np.arange(4)*(nx-1)/4
    #yc=xc*0+cel[i,1]
    #slat=ylat[yc]
    #slon=xlon[xc]

    #;  Convert slat slon to elat elon with given new pole
    #   elat, elon, sxc, syc = steromap( slat, slon, plat, plon, Onecl=True )

    #   if( (elat[0] >= 0.0) and (rngsxy[0] < sxc[0] < rngsxy[1])
    #                        and (rngsxy[2] < syc[0] < rngsxy[3]) ):
    #       nvrts.append(list(zip(sxc,syc)))
    #       ncels.append(i)

    #   if( (elat[0] <  0.0) and (pltype == 'Global') ):
    #       svrts.append(list(zip(sxc,syc)))
    #       scels.append(i)

    #  Set plot size and limits and message out anchor point.
    rdpols=[radius, pangle, plon, plat]
    config=np.array([rdpols, sztpxy, rngsxy, clrbxy, ngabjm])
    pzfile=DatGMC+'/' + ModlName + '_Vrts'+pltype[0:4]+'.npz'
    #pzfile=DatGMC+'/' + ModlName + '_Vrts.npz'

    #  Store selected north and south verts and cell numbers for swh plots.
    #  Use the np.savez to save 3/5 variables in one file.  JGLi22Feb2019 
    if( pltype == 'Global' ):
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
                          svrt=svrts, scel=scels)
    else:
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 
    #  These variables can be loaded back using
    #   vrtcls = np.load(DatGMC+'S625Vrts'+pltype[0:4]+'.npz')
    #   nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']; config=vrtcls['cnfg']; 
    #   svrts = vrtcls['svrt'] ; scels = vrtcls['scel']

    #  Draw your selected grid plot.
    psfile=Wrkdir + '/' + ModlName + '_' + pltype[0:4] + 'grd' + figtype 

    if( pltype == 'Global'):
        smcglobl( cel, nvrts,ncels,svrts,scels,colrs, config,
                 mdlname= ModlName, Arctic=Arctic, buoys=None, psfile=psfile)

    else:
        smclocal( cel, nvrts,ncels,colrs,config, Arctic=Arctic, 
                 mdlname= ModlName, buoys=None, psfile=psfile,
                 paprorn=papror)

    print( "[INFO] Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

    #End of pltSMCgrid program
    return


#---
# main program

ModlName = sys.argv[1]
Wrkdir   = sys.argv[2]

# read the grid metadata from the grid namelist file
gridmeta = readGridnml(Wrkdir,fname='smcGrid.nml')

# plot the grid
pltSMCgrid(ModlName, gridmeta['nlevs'], 
           gridmeta['dx'], gridmeta['dy'], 
           gridmeta['x0'], gridmeta['y0'], 
           gridmeta['nx'], gridmeta['ny'], 
           Wrkdir, Cel_file=gridmeta['grdfile'],
           Arc_file=gridmeta['arcfile'], figtype='.ps')
