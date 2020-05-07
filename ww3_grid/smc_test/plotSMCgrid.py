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
# plotSMCgrid.py
#
# PURPOSE:
#  Functions library for generating spherical multiple-cell (SMC) grid plots
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
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import pandas as pd

"""
  Included Functions:

    readcell(celfiles)

    readtext(textfile, skiprows=[0])

    rgbcolor( colrfile )

    ll2eqdeg( SLat, SLon, Polat, Polon )

    steromap( SLat, SLon, Polat, Polon, Pangl=90.0, Onecl=False )

    scale_depth( ncstr=0, nclrm=136 )

    colrboxy(cboxy, colrmap, marks, ncstr=0, nclrm=256)

    smcglobl(cel, nvrts,ncels,svrts,scels,colrs, config, Arctic=True, 
             mdlname='SMC', buoys=None, psfile='output.ps', paprtyp='a3')

    smclocal(cel, verts, ncels, colrs, config, Arctic=False, 
             mdlname='SMC', buoys=None, psfile='output.ps', 
             paprorn='portrait', paprtyp='a3')

    smcglobl(cel, nvrts,ncels,svrts,scels,colrs, config, Arctic=True, 
             mdlname='SMC', buoys=None, psfile='output.ps', paprtyp='a3')

    swhlocal(swhs,verts,ncels,colrs,config,
        mdlname='SMC',datx='2018',psfile='output.ps', 
        paprorn='portrait', paprtyp='a3')
"""

def readcell(celfiles):
    """Read SMC grid cell arrays from text files.

        Read SMC grid cell files and return count line as headrs
        Function readcell(filestring) read the cell file specified by
        the filestring, which includes path/file_name.

        If Arctic cell file is also needed, the filestring will be a
        list of both global and Arctic cell files.  It will merge the 
        two parts together

        JGLi18Feb2019 
        Modified to use numpy.genfromtxt.   JGLi10Feb2020"""

    ##  Input celfiles as a list even if there is only one file.
    ##  For instance celfiles=['path/celfile.dat', 'path/arcfile.dat']
    nfls = 0 
    for celfile in celfiles:
        print( " Read cel from ", celfile)
        archd=open(celfile, 'r') 
        hdlin=archd.readline()
        archd.close()
         
        celin=np.genfromtxt(celfile, skip_header=1)
        nfls += 1

        if( nfls <= 1):
            headrs = [ hdlin ]
            cel = celin.copy()
        else:
            hdlin = hdlin + ' %d' %np.min(celin[:,1]) # add boundary cell j-value
            headrs.append( hdlin )
            cel = np.vstack( (cel, celin) )
        
    return headrs, cel.astype(np.int)


def readtext(textfile, skiprows=[0]):
    """
    Read data arrays from a tabled text file. Missing values in the
     table will be filled as none.

     Function readtext(filestring) read the cell file specified by
     the filestring, which includes path/file_name.txt

              JGLi18Feb2019"""

    ##  Assume to skip the first line [0]. If more than one
    ##  lines to skip, specifiy the line indexes in the
    ##  parameter list skiprows=[0,1], for first two lines.
    ##  But only the last skipped line is returned as a list.

    print (" Read table from ", textfile)
    archd=open(textfile, 'r') 
    #  Assume last line in skiprows contains the header parameters.
    for i in range(len(skiprows)):
        hdrskp=archd.readline()

    #  Split the last one of skipped lines as header parameters.
    hdlist = hdrskp.split()
    archd.close()

    arcel=pd.read_csv(textfile, sep='\s+',skiprows=skiprows,header=None)
    table=arcel.values

    return hdlist, table


def rgbcolor( colrfile ):
    """
    #;;  rgbcolor function to set up the rainbow color with 0 to 255 index. 
    #;;  First Created on  18 Dec 2018  by Jian-Guo Li
    #;;  Last Modified on  20 Dec 2018  by Jian-Guo Li
    #
    # name:     rgbcolor
    #
    # purpose:  read rgbspectrum.dat color table file and setup as a colormap in Python.
    #
    # usage:    colrs = rgbcolor() 
    #
    # input:    rgbspectrum.dat  color table file in the current directiory.
    #           readcell.py  
    #
    # output:   colrs colormap object of 256 tuples, each may be used as colrs(i) 
    #
    """

    nclrs, clsrgb = readtext( colrfile )
    clsnm = clsrgb/float( clsrgb.max() ) 
    colrs = mpl.colors.ListedColormap(clsnm, name='Spectra')

    nclrm = colrs.N  ## colrs.N returns the total number of colors 256.
    normd = mpl.colors.Normalize(vmin=0, vmax=nclrm-1)

    return colrs


def ll2eqdeg( SLat, SLon, Polat, Polon ):
    """  ll2eqdeg procedure converts standard latitude longitude to rotated 
    #    or equatorial latitude longitude with a given new pole position.
    #    First created on  4 Dec 2018  by Jian-Guo Li
    #    Last modified on 10 Dec 2018  by Jian-Guo Li
    #
    # name:     ll2eqdeg
    #
    # purpose:  Converts standard lat/lon to rotated lat/lon as NumPy arrays.
    #
    # usage:    ELat, ELon = ll2eqdeg( SLat, SLon, PoLat, PoLon )
    #
    # input:    SLat,  SLon  --- Standard lat lon as NumPy array in deg
    #           Polat, Polon --- New North Pole standard lat/lon in deg, scalors.
    # output:   ELat,  ELon  --- Corresponding lat lon in rotated grid in deg, NumPy array
    #
    """

    #  Check Input SLat, Slon elements, should be equal
    nlat=len(SLat)
    nlon=len(SLon)

    if nlat != nlon: 
        print( ' SLat and SLon elements should have equal elements!' )
        print( ' SLat and SLon elements are', nlat, nlat )
        return 

    #  No need to calculate if North Pole unchanged
    if (Polat == 90.0 and Polon == 0.0): 
        return (SLat, SLon)

    #  Constants 
    D2Rad=np.pi/180.0
    R2Deg=180.0/np.pi

    #  Make Pole longitude within range -180 to 180
    PoLon_orig=Polon
    if Polon > 180.0:  Polon=Polon - 360.0

    #  Sine and cosine of PoLat
    if Polat >= 0.0: 
        sinpolat= np.sin(Polat*D2Rad)
        cospolat= np.cos(Polat*D2Rad)
    else: 
        sinpolat=-np.sin(Polat*D2Rad)
        cospolat=-np.cos(Polat*D2Rad)

    #  Shifting SLat to opposite rotated pole side. 
    ZeroLon = Polon + 180.0
    ALon = SLon - ZeroLon

    #  Converting ALon range to be -180 to 180
    lonind= ALon >= 180.0 
    if( lonind.any() ): ALon[lonind]=ALon[lonind]-360.0
    lonind= ALon < -180.0 
    if( lonind.any() ): ALon[lonind]=ALon[lonind]+360.0

    #  Conversion of SLat/SLon to ELat/ELon
    Apt_Ang = - cospolat*np.cos(SLat*D2Rad)*np.cos(ALon*D2Rad) \
              + sinpolat*np.sin(SLat*D2Rad)

    Apt_Ang[ Apt_Ang >  1.0 ]= 1.0 
    Apt_Ang[ Apt_Ang < -1.0 ]=-1.0 

    ELatRad = np.arcsin(Apt_Ang)
    cosELat = np.cos(ELatRad)
    cosELon = sinpolat*np.cos(SLat*D2Rad)*np.cos(ALon*D2Rad) \
             +cospolat*np.sin(SLat*D2Rad)

    ELat = ELatRad*R2Deg
    ELon = np.zeros(nlon)
    Tmpr = np.zeros(nlon)

    #  Only set Elon where cosELat is non zero
    latind = cosELat > 0.0 
    if( latind.any() ): 
        Tmprat=cosELon[latind]/cosELat[latind]

        Tmprat[ Tmprat >  1.0 ] = 1.0
        Tmprat[ Tmprat < -1.0 ] =-1.0

        ELon[latind] = R2Deg*np.arccos(Tmprat)

    #  Change ELon sign to negative if ALon is negative. 
    lonind= ALon < 0.0 
    if( lonind.any() ): ELon[lonind]= - ELon[lonind]

    #print( '... Finishing ll2eqdeg ...') 
    return (ELat, ELon)


def steromap( SLat, SLon, Polat, Polon, Pangl=90.0, Onecl=False ):
    """
    #;;  steromap procedure converts standard latitude longitude to rotated 
    #;;  or equatorial latitude longitude with a given new pole position and
    #;;  map the rotated lat/lon with a sterographic projection.
    #;;  Created on  14 May 2014  by Jian-Guo Li
    #;;  Modified on 11 Nov 2014  by Jian-Guo Li
    #;;  Converted into a Python function on 5 Dec 2018 by Jian-Guo Li
    #
    # name:     steromap
    #
    # purpose:  Converts standard lat/lon to rotated lat/lon and sterographic map
    #
    # usage:    ELat, ELon, SXxc, SYyc = steromap( SLat, SLon, Polat, Polon, Pangl=Pangl, Onecl=True )
    #
    # input:    SLat, SLon --- Standard lat lon as ndarray in deg
    #           Polat, Polon --- New North Pole position in standard lat/lon in deg, scalors.
    #           Pangl --- scalor angle in deg from rotated Pole so its projected radius is 10 unit.
    #           Onecl --- if True, keep all points in one hemisphere as one cell. 
    # output:   ELat, ELon --- Corresponding lat lon in rotated grid in deg, as ndarray.
    #           SXxc, SYyc --- Corresponding projected x/y coordinates in range [-10, 10], as ndarray.
    #
    """

    # PRO  STEROMAP, SLat, SLon, ELat, ELon, SXxc, SYyc, Polat=Polat, Polon=Polon, Pangl=Pangl, Onecl=Onecl
 
    #  Check Input SLat, Slon elements, should be equal
    nlat=len(SLat)
    nlon=len(SLon)

    if nlat != nlon: 
        print( ' SLat and SLon elements should have equal elements!' )
        print( ' SLat and SLon elements are', nlat, nlat )
        return 

    if Pangl <= 0.0: 
        print( ' Pangl has to be 0 < Pangl <= 90.0!' )
        print( ' Input Pangl is equal to', Pangl )
        return 

    #  No need to calculate if North Pole unchanged
    if (Polat == 90.0 and Polon == 0.0): 
        ELat = SLat
        ELon = SLon
    else:
    #  Convert slat slon to elat elon with given new pole
        ELat, ELon = ll2eqdeg(SLat, SLon, Polat, Polon)

    #  Projection parameters
    #  Constants 
    d2rad=np.pi/180.0
    r2deg=180.0/np.pi
    radius=10.0

    #  Number of radius of projection distance
    pnrds=4.0

    #  Adjusted projected radius so that projected edge radius 
    #  to be equal to the original radius. 
    PRadus=radius*(pnrds - 1.0 + np.cos(Pangl*d2rad))/np.sin(Pangl*d2rad)

    #  Generate projecting coordiantes. ndarray operations.
    pradmp=PRadus*np.cos(ELat*d2rad)/(pnrds - 1.0 + np.sin(np.absolute(ELat)*d2rad))
    SYyc = pradmp*np.cos(ELon*d2rad)
    SXxc =-pradmp*np.sin(ELon*d2rad)

    #  Check if any point is on the southern hemisphere
    indx = ELat < 0.0 
    if( indx.any() ):     ## Any element of this boolean array is true
        if( Onecl ):
            #  If it is for one single cell projection, keep all in one hemisphere
            SXxc = -SXxc 
        else:
            #  Reverse southern hemisphere x-coordinate for individual points.
            SXxc=np.where(indx, -SXxc, SXxc)

    #print('... Finishing steromap.np ...')
    return ( ELat, ELon, SXxc, SYyc )


def scale_depth( ncstr=0, nclrm=136 ):
    """
    #;;  scale_depth function to set up the color scale for depth plot, using the
    #;;  first 130 colors from a given color table of 256 colors.
    #;;  First Created on  22 Feb 2019  by Jian-Guo Li
    #;;  Last Modified on  26 Feb 2019  by Jian-Guo Li
    #
    # usage:    depth, factr, cstar, marks, ckeys = scale_depth( nclrm=131 )
    #
    # Bathymetry depth (positive) will need to be converted by
    #    ndeps=np.rint( (cstar - np.log10(depth))*factr ).astype(np.int)
    # to be consistent with the color key. 
    #
    """

    depth=np.array([10000,1000,100,10,0,-10], dtype=np.int)
    cstar=np.log10(depth[0]+1000.0)
    factr=(nclrm - 1)/cstar
    marks= ncstr + np.rint( (cstar - np.log10(depth+11))*factr ).astype(np.int)
    #print( nclrm, ' colors with depth marks at ', depth[::-1])

    return ( depth, factr, cstar, marks, ncstr, nclrm )


def scale_swh( ncstr=0, nclrm=256 ):
    """
    #;;  scale_swh function to set up the color scale for SWH plot, using given
    #;;  color table with 0 to 255 index. 
    #;;  First Created on  22 Feb 2019  by Jian-Guo Li
    #;;  Last Modified on  26 Feb 2019  by Jian-Guo Li
    #
    # usage:    waveht, factor, residu, marks, ckeys = scale_swh( nclrm=colrs.N )
    #
    # For consistency wave height should be converted by
    #   nswh=np.rint( factor*np.log(waveht+residu) ).astype(np.int)
    #
    """

    waveht=np.array([0,1,2,4,8,16,32], dtype=np.int)
    factor=(nclrm - 2)/np.log(35.0)
    residu=np.exp(5.0/factor)
    marks=ncstr + np.rint( factor*np.log(waveht+residu) ).astype(np.int)
    print( ncstr, nclrm, ' colors with waveheight marks at ', waveht[:] )

    return ( waveht, factor, residu, marks, ncstr, nclrm )


def colrboxy(cboxy, colrmap, marks, ncstr=0, nclrm=256):
    """
    Generate a colour bar polycollection with given key box 
    cboxy=[x, y, dx, dy] and colormap, plus marks.
    The colour index arrays will generated out of given colormap.

                  JGLi01Mar2019 
    """

    ##  Input variables include x, y as Numpy Arrays or lists for the 
    ##  colorbar locations, the colormap to be plotted, integer list 
    ##  or Numpy Array marks for key scale ticks as enlarged polygon. 
    ##  Tick labels should be added outside this program when the 
    ##  colorbar PolyCollection is drawn in a corresponding plot. 

    ## Workout color bar orientaion by x and y sizes
    ## and generate poly verts and colors. 
    x0=cboxy[0]
    y0=cboxy[1]
    bx=cboxy[2] 
    by=cboxy[3] 
    verts = []
    pcolr = []
    ic = ncstr

    if( bx > by ):  ## Horizontal color bar
        dx = bx/nclrm
        xkeys=np.arange(nclrm)*dx + x0
        ykeys=np.array([y0, y0+by])
        syc=[y0,y0,y0+by,y0+by]
        for xi in xkeys:
            sxc=[xi, xi+dx, xi+dx, xi]
            verts.append(list(zip(sxc,syc)))
            pcolr.append(colrmap(ic))
            ic += 1

        dm = 1.1*by
        sym=[y0,y0,y0+dm,y0+dm]
        for i in marks:
            xi=xkeys[i]
            sxc=[xi, xi+dx, xi+dx, xi]
            verts.append(list(zip(sxc,sym)))
            pcolr.append(colrmap(i))
 
    else:           ## Vertical color bar
        dy = by/nclrm 
        xkeys=np.array([x0, x0+bx])
        ykeys=np.arange(nclrm)*dy + y0
        sxc=[x0,x0+bx,x0+bx,x0]
        for yj in ykeys:
            syc=[yj, yj, yj+dy, yj+dy]
            verts.append(list(zip(sxc,syc)))
            pcolr.append(colrmap(ic))
            ic += 1

        dm = 1.1*bx
        sxm=[x0,x0+dm,x0+dm,x0]
        for j in marks:
            yj=ykeys[j]
            syc=[yj, yj, yj+dy, yj+dy]
            verts.append(list(zip(sxm,syc)))
            pcolr.append(colrmap(j))

    ## Generate PolyCollection
    cbarpoly = PolyCollection(verts)
    cbarpoly.set_color(pcolr)

    return ( xkeys, ykeys, cbarpoly )


def smcglobl(cel, nvrts, ncels, svrts, scels, colrs, config, Arctic=None, 
             mdlname='SMC', buoys=None, psfile='output.ps', paprtyp='a3'):
    """
    ##  The smcglobl function plots a global view of a given smc grid.
    ##  cel is the full cell array in shape(nc, 5).
    ##                                JGLi04Mar2019
    """

    ##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

    ##  Global plot configuration parameters.
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]
    if( Arctic is not None ): ncabgm=config[4]

    ##  Maximum mapping radius.
    radius=rdpols[0]
    pangle=rdpols[1]
    plon  =rdpols[2]
    plat  =rdpols[3]

    ##  Outline circle at equator
    ciran=np.arange(1081)*d2rad/3.0
    xcirc=radius*np.cos(ciran)
    ycirc=radius*np.sin(ciran)

    ##  Define depth to color index conversion parameters and marks.
    ##  Use only first 131 colors in colrs(0:255). 
    depth, factr, cstar, marks, ncstr, nclrm = scale_depth(nclrm=131)

    ##  Cell color is decided by its depth value, dry cells use default 0 color.
    nc = cel.shape[0]
    ndeps = np.zeros( (nc), dtype=np.int )
    for j in range(nc):
        if( cel[j,4] > 0 ):
            ndeps[j] = ncstr + np.rint( (cstar-np.log10(cel[j,4]))*factr ).astype(np.int)
    #ndeps = ncstr + np.rint( (cstar-np.log10(cel[:,4]))*factr ).astype(np.int)

    if( Arctic is not None ):
        na=int(ncabgm[1])
        nb=int(ncabgm[2])
        jm=int(ncabgm[3])

    ##  Set up first subplot and axis for northern hemisphere
    print (" Drawing north hemisphere cells ... ")
    fig=plt.figure(figsize=sztpxy[0:2])
    ax1=fig.add_subplot(1,2,1)
    ax1.set_aspect('equal')
    ax1.set_autoscale_on(False)
    ax1.set_axis_off()
    plt.subplots_adjust(left=0.02,bottom=0.02,right=0.98,top=0.95)
    plt.subplots_adjust(wspace=0.02, hspace=0.01)

    xprop={'xlim':rngsxy[0:2], 'xlabel':'Nothing'}
    yprop={'ylim':rngsxy[2:4], 'ylabel':'Nothing'}
    ax1.set(**xprop)
    ax1.set(**yprop)

    ax1.plot(xcirc, ycirc, 'b-')

    ## Select cells for one subplot and create verts and pcolr. 
    pcface = []
    pcedge = []
    for i in ncels:
        if( Arctic is not None ) and ( cel[i,1] == Arctic ):
            # Mark the arctic map-east reference region by first ring of its boundary cells
            pcface.append(colrs(246))
            pcedge.append(colrs(246))
        #elif( Arctic is not None ) and ( cel[i,1] == jm ):
        #    # Mark the arctic cells
        #    pcface.append(colrs(168))
        #    pcedge.append(colrs(168))
        else:
            pcface.append(colrs(255))
            pcedge.append(colrs(ndeps[i]))
                
    ##  Draw this hemisphere cells
    smcpoly = PolyCollection(nvrts)
    smcpoly.set_facecolor(pcface)
    smcpoly.set_edgecolor(pcedge)
    smcpoly.set_linewidth( 0.2 )
    ax1.add_collection(smcpoly)

    ##  Draw colorbar for ax1.
    xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks, nclrm=nclrm)
    ax1.add_collection(clrply)
    for i in range(len(depth)):
        m = marks[i]
        plt.text(xkeys[m], ykeys[0]+1.15*(ykeys[1]-ykeys[0]), str(depth[i]),
            horizontalalignment='center', fontsize=11, color='b' )
            #rotation=-90,verticalalignment='center', fontsize=11, color='b' )

    plt.text(xkeys[marks[0]], ykeys[0]+2.2*(ykeys[1]-ykeys[0]), 'Depth m',
            horizontalalignment='left', fontsize=15, color='k' )
           #rotation=-90,verticalalignment='center', fontsize=15, color='k' )

    # Overlay buoy sites on grid map if buoy file is provided.
    if(buoys is not None):
        hdr, buoyll = readtext(buoys)
        nmbu=long(hdr[0])
        buoyids=buoyll[:,0].astype(str)
        buoylat=buoyll[:,1].astype(np.float)
        buoylon=buoyll[:,2].astype(np.float)

        #; Convert slat slon to elat elon with given new pole
        elat,elon,sxc,syc = steromap(buoylat,buoylon,plat,plon,Pangl=pangle)

        #; Mark buoy position on map
        print (' Selected buoys on north hemisphere ...')
        for i in range(nmbu):
            if( (elat[i] >= 0.0) and (rngsxy[0] < sxc[i] < rngsxy[1])
                                 and (rngsxy[2] < syc[i] < rngsxy[3]) ):
                print (' {:6} {:8.3f} {:8.3f}'.format( buoyids[i], buoylat[i], buoylon[i] ) )
                txtsz=int( abs(np.sin(elat[i]*d2rad)*12.0) )
                plt.text(sxc[i], syc[i], 'r',  fontsize=txtsz,
                     horizontalalignment='center', color='r' )
                plt.text(sxc[i], syc[i], '.',  fontsize=txtsz*2,
                     horizontalalignment='center', color='r' )


    ##  Southern hemisphere
    print (" Drawing south hemisphere cells ... ")
    ax2=fig.add_subplot(1,2,2)
    ax2.set_aspect('equal')
    ax2.axis('off') 
    plt.subplots_adjust(left=0.02,bottom=0.02,right=0.98,top=0.95)
    plt.subplots_adjust(wspace=0.02, hspace=0.01)
    ax2.set(**xprop)
    ax2.set(**yprop)

    ax2.plot(xcirc, ycirc, 'b-')

    ## Select cells for one subplot and create verts and pcolr. 
    pcface = []
    pcedge = []
    for i in scels:
        pcface.append(colrs(255))
        pcedge.append(colrs(ndeps[i]))

    ## Generate polygon collections for southern hemisphere
    smcpoly = PolyCollection(svrts)
    smcpoly.set_facecolor(pcface)
    smcpoly.set_edgecolor(pcedge)
    smcpoly.set_linewidth( 0.2 )

    ax2.add_collection(smcpoly)
  
    ##  Put cell information inside plot
    tpx=sztpxy[2] 
    tpy=sztpxy[3]
    dpy= 0.6 

    plt.text(tpx, tpy+dpy*0.5, mdlname+' Grid',  
             horizontalalignment='center', fontsize=15, color='k' )
    plt.text(tpx, tpy+dpy*2.0, 'NC='+str(nc), 
             horizontalalignment='center', fontsize=13, color='r' )
    if( Arctic is not None ):
        plt.text(tpx, tpy+dpy*3.0, 'NA='+str(na), 
             horizontalalignment='center', fontsize=13, color='r' )
        plt.text(tpx, tpy+dpy*4.0, 'NB='+str(nb), 
             horizontalalignment='center', fontsize=13, color='r' )

    ##  Draw colorbar for ax2.
    xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks, nclrm=nclrm)
    ax2.add_collection(clrply)
    for i in range(len(depth)):
        m = marks[i]
        plt.text(xkeys[m], ykeys[0]+1.18*(ykeys[1]-ykeys[0]), str(depth[i]),
            horizontalalignment='center', fontsize=11, color='b' )
            #verticalalignment='center', fontsize=11, color='b' )
            #rotation=-90,verticalalignment='center', fontsize=11, color='b' )

    plt.text(xkeys[marks[-1]], ykeys[0]+2.2*(ykeys[1]-ykeys[0]), 'Depth m',
            horizontalalignment='right', fontsize=15, color='k' )
            #rotation=-90,verticalalignment='center', fontsize=15, color='k' )

    # Overlay buoy sites on grid map if buoy file is provided.
    if(buoys is not None):

    #; Mark buoy position on map
        print (' Selected buoys on south hemisphere ...')
        for i in range(nmbu):
            if( (elat[i] < 0.0) and (rngsxy[0] <-sxc[i] < rngsxy[1])
                                and (rngsxy[2] < syc[i] < rngsxy[3]) ):
                print (' {:6} {:8.3f} {:8.3f}'.format( buoyids[i], buoylat[i], buoylon[i] ) )
                txtsz=int( abs(np.sin(elat[i]*d2rad)*12.0) )
                plt.text(-sxc[i], syc[i], 'r',  fontsize=txtsz,
                     horizontalalignment='center', color='r' )
                plt.text(-sxc[i], syc[i], '.',  fontsize=txtsz*2,
                     horizontalalignment='center', color='r' )

    ##  Refresh subplots and save them.
    plt.subplots_adjust(wspace=0.02, hspace=0.01)

    plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
                orientation='landscape',papertype=paprtyp,format='ps')

    # End of smcglobl plot function


def smclocal(cel, verts, ncels, colrs, config, Arctic=None, 
             mdlname='SMC', buoys=None, psfile='output.ps', 
             paprorn='portrait', paprtyp='a3'):
    """
    ##  The smclocal function plots a local region of a given smc grid.
    ##  cel is the full cell array in shape(nc, 5).
    ##                                JGLi04Mar2019
    """

    #  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

    #  Local plot configuration parameters
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]
    if( Arctic ): ncabgm=config[4]

    # Maximum mapping radius.
    radius=rdpols[0]
    pangle=rdpols[1]
    plon  =rdpols[2]
    plat  =rdpols[3]

    # Define depth to color index conversion parameters and marks.
    # Use only first 131 colors in colrs(0:255). 
    depth, factr, cstar, marks, ncstr, nclrm = scale_depth(nclrm=136)

    # Cell color is decided by its depth value.
    nc = cel.shape[0]
    ndeps = np.zeros( (nc), dtype=np.int )
    for j in range(nc):
        if( cel[j,4] > -11 ): 
            ndeps[j] = ncstr + np.rint( (cstar-np.log10(cel[j,4]+11))*factr ).astype(np.int)
    #ndeps = ncstr + np.rint( (cstar-np.log10(cel[:,4]))*factr ).astype(np.int)

    if( Arctic is not None ):
        na=int(ncabgm[1])
        nb=int(ncabgm[2])
        jm=int(ncabgm[3])

    # Workout output file format from its extension
    #for i in np.arange(len(psfile))[::-1]:
    #    if( psfile[i] == '.' ):
    #        psfmt = psfile[i+1:]
    #        break
    #print " Output file format will be "+psfmt
    # Python plt.savefig will do this automatically.

    # Use selected cells to draw the plot.
    fig=plt.figure(figsize=sztpxy[0:2])
    ax=fig.add_subplot(1,1,1)
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    ax.set(**xprop)
    ax.set(**yprop)
    ax.set_aspect('equal')
    ax.set_autoscale_on(False)
    ax.set_axis_off()
    plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

    # Create color array for this plot
    pcface = []
    pcedge = []
    for i in ncels:
        if( Arctic is not None ) and ( cel[i,1] == Arctic ):
            # Mark the arctic map-east reference region by first ring of its boundary cells
            pcface.append(colrs(246))
            pcedge.append(colrs(246))
        #elif( Arctic is not None ) and ( cel[i,1] == jm ):
        #    # Mark the arctic cells
        #    pcface.append(colrs(168))
        #    pcedge.append(colrs(168))
        else:
            pcedge.append(colrs(ndeps[i]+10))
            pcface.append(colrs(ndeps[i]))
            #pcface.append(colrs(255))
            #pcedge.append(colrs(ndeps[i]))

    # Create PolyCollection from selected verts and define edge and face color.
    polynorth = PolyCollection(verts)
    polynorth.set_facecolor(pcface)
    polynorth.set_edgecolor(pcedge)
    polynorth.set_linewidth( 0.2 )

    # Draw the selected cells as colored polygons.
    ax.add_collection(polynorth)

    # Draw colorbar inside plot.
    xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks, nclrm=nclrm)
    ax.add_collection(clrply)
    dkx=clrbxy[2]; dky=clrbxy[3]
    for i in range(len(depth)):
        m = marks[i]
        if( dkx < dky ):    
            plt.text(xkeys[0]+1.15*dkx, ykeys[m], str(depth[i]),
                verticalalignment='center', fontsize=11, color='b' )
        else:
            plt.text(xkeys[m], ykeys[0]+1.15*dky, str(depth[i]),
                horizontalalignment='center', fontsize=11, color='b' )

    if( dkx < dky ):    
        plt.text(xkeys[0]+1.9*dkx, ykeys[marks[2]], 'Depth m',
                 rotation=-90,verticalalignment='center', fontsize=15, color='k' )
    else:
        plt.text(xkeys[marks[2]], ykeys[0]+1.9*dky, 'Depth m',
                 rotation=0,horizontalalignment='center', fontsize=15, color='k' )

    # Put cell information inside plot
    tpx=sztpxy[2] 
    tpy=sztpxy[3]
    dpy=-0.6 

    plt.text(tpx, tpy+dpy*1, mdlname+' Grid',  
             horizontalalignment='center', fontsize=15, color='k' )
    plt.text(tpx, tpy+dpy*2, 'NC='+str(nc), 
             horizontalalignment='center', fontsize=13, color='r' )
    if( Arctic is not None ):
        plt.text(tpx, tpy+dpy*3, 'NA='+str(na), 
             horizontalalignment='center', fontsize=13, color='r' )
        plt.text(tpx, tpy+dpy*4, 'NB='+str(nb), 
             horizontalalignment='center', fontsize=13, color='r' )

    # Overlay buoy sits on grid map if buoy file is provided.
    if(buoys is not None):
        hdr, buoyll = readtext(buoys)
        nmbu=int(hdr[0])
        buoyids=buoyll[:,0].astype(str)
        buoylat=buoyll[:,1].astype(np.float)
        buoylon=buoyll[:,2].astype(np.float)

        # Convert slat slon to elat elon with given new pole
        elat,elon,sxc,syc = steromap(buoylat,buoylon,plat,plon,Pangl=pangle)

        # Mark buoy position on map
        print (' Selected buoys in this plot:')
        for i in range(nmbu):
            if( (elat[i] >= 25.0) and (rngsxy[0] < sxc[i] < rngsxy[1])
                              and (rngsxy[2] < syc[i] < rngsxy[3]) ):
                print (' {:6} {:8.3f} {:8.3f}'.format( buoyids[i], buoylat[i], buoylon[i] ))
                txtsz=int( abs(np.sin(elat[i]*d2rad)*12.0) )
                plt.text(sxc[i], syc[i], 'r',  fontsize=txtsz,
                     horizontalalignment='center', color='r' )
                plt.text(sxc[i], syc[i], '.',  fontsize=txtsz*2,
                     horizontalalignment='center', color='r' )

    # Save plot as ps file
    print (" Save the smc grid local plot ... " )
    plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
                orientation=paprorn,papertype=paprtyp,format='ps')

    #  End of smclocal plot function. ##


"""
##  Draw global swh plot with given polycollections and cell
##  array. Output as A3 ps file.        JGLi28Feb2019
##
"""

def swhglobl(swhs, nvrts, ncels, svrts, scels, colrs, config,
             mdlname='SMC', datx='2018010106',
             psfile='output.ps', paprtyp='a3'):
    """
    ##  Draw global swh plot with given polycollections and cell
    ##  array. Output as A3 ps file.        JGLi28Feb2019
    ##
    """

    # Degree to radian conversion parameter.
    d2rad=np.pi/180.0

    # Global plot configuration parameters.
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]
    ncabgm=config[4]

    # Maximum mapping radius.
    radius=rdpols[0]
    pangle=rdpols[1]
    plon  =rdpols[2]
    plat  =rdpols[3]

    # Outline circle at equator
    ciran=np.arange(1081)*d2rad/3.0
    xcirc=radius*np.cos(ciran)
    ycirc=radius*np.sin(ciran)

    # Set up wave height scale and marks with colrs.N.
    # colrs.N returns colrs' total number of colors 256.
    waveht, factor, residu, marks, ncstr, nclrm = scale_swh(nclrm=colrs.N)

    resmn1=residu - 1.0
    nswh0=ncstr+ int( factor*np.log(-resmn1 + residu) )
    #print (' factor, residu, resmn1, nswh0 = {} {} {} {: d}' 
    #.format(factor, residu, resmn1, nswh0) ) 

    # Some constant variables for plots.
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    
    if True:
        # Work out max and min values, excluding missing data (-999.0)
        cmax = swhs.max()
        cmin = swhs[ swhs > -999.0 ].min()
        print ( ' swh range %f, %f' % (cmin, cmax) )
        cmxs = 'SWHmx = %6.2f m' % cmax
        cmns = 'SWHmn = %10.3E' % cmin

        # Reset missing values (-999.0) to be -resmn1 
        swhs[ swhs < -resmn1] = -resmn1

        # Trim large values into plot range if any
        swhs[ swhs > 32.0 ] = 32.0 

        # Convert swhs with logarithm scale.
        icnf = np.rint( factor*np.log(swhs+residu) )
        nswh = np.array( icnf, dtype=np.int )

        print (" Drawing "+psfile)
        # Set up first subplot and axis for northern hemisphere
        fig=plt.figure(figsize=sztpxy[0:2])
        ax1=fig.add_subplot(1,2,1)
        ax1.set_aspect('equal')
        ax1.set_autoscale_on(False)
        ax1.set_axis_off()
        plt.subplots_adjust(left=0.02,bottom=0.02,right=0.98,top=0.95)
        plt.subplots_adjust(wspace=0.02, hspace=0.01)

        ax1.set(**xprop)
        ax1.set(**yprop)

        ax1.plot(xcirc, ycirc, 'b-')

        # Use loaded verts to setup PolyCollection.
        polynorth = PolyCollection(nvrts)

        # Create color array for this plot
        pcface = []
        pcedge = []
        for i in ncels:
            if( nswh[i] == nswh0 ): 
                pcface.append(colrs(255)) 
                pcedge.append(colrs(0)) 
            else:
                pcface.append(colrs(nswh[i]))
                pcedge.append(colrs(nswh[i]))
               

        #smcpoly.set_color(pcolr)    ## This line defines both edge and face color.
        polynorth.set_facecolor(pcface)
        polynorth.set_edgecolor(pcedge)
        polynorth.set_linewidth( 0.2 )
        ax1.add_collection(polynorth)  

        # Draw colorbar for ax1.
        xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks)
        ax1.add_collection(clrply)
        for i in range(len(waveht)):
            m = marks[i]
            plt.text(xkeys[m], ykeys[0]+1.18*(ykeys[1]-ykeys[0]), str(waveht[i]),
            horizontalalignment='center', fontsize=11, color='b' )
            #rotation=-90,verticalalignment='center', fontsize=11, color='b' )

        plt.text(xkeys[marks[0]], ykeys[0]+2.2*(ykeys[1]-ykeys[0]), 'SWH m',
            horizontalalignment='left', fontsize=15, color='k' )
            #rotation=-90,verticalalignment='center', fontsize=15, color='k' )


        # Southern hemisphere subplot.
        ax2=fig.add_subplot(1,2,2)
        ax2.set_aspect('equal')
        ax2.axis('off')
        plt.subplots_adjust(left=0.02,bottom=0.02,right=0.98,top=0.95)
        plt.subplots_adjust(wspace=0.02, hspace=0.01)
        ax2.set(**xprop)
        ax2.set(**yprop)

        ax2.plot(xcirc, ycirc, 'b-')

        # Use loaded verts to set up PolyCollection. 
        polysouth = PolyCollection(svrts)

        # Create color array for this plot
        pcface = []
        pcedge = []
        for i in scels:
            if( nswh[i] == nswh0 ):
                pcface.append(colrs(255))
                pcedge.append(colrs(0))
            else:
                pcface.append(colrs(nswh[i]))
                pcedge.append(colrs(nswh[i]))
   
        #   smcpoly.set_color(pcolr)    ## This line defines both edge and face color.
        polysouth.set_facecolor(pcface)
        polysouth.set_edgecolor(pcedge)
        polysouth.set_linewidth( 0.2 )
        ax2.add_collection(polysouth)


        # Put statistic information inside subplot ax2
        tpx=sztpxy[2] 
        tpy=sztpxy[3]
        dpy= 0.6 

        plt.text(tpx, 9.0, mdlname+' SWH',  
             horizontalalignment='center', fontsize=19, color='r' )
        plt.text(tpx, tpy+dpy*1.0, cmns,
             horizontalalignment='center', fontsize=15, color='b' )
        plt.text(tpx, tpy+dpy*2.0, cmxs, 
             horizontalalignment='center', fontsize=15, color='r' )
        plt.text(tpx, tpy+dpy*3.0, datx,
             horizontalalignment='center', fontsize=17, color='k' )

        # Draw colorbar for ax2.
        xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks)
        ax2.add_collection(clrply)
        for i in range(len(waveht)):
            m = marks[i]
            plt.text(xkeys[m], ykeys[0]+1.18*(ykeys[1]-ykeys[0]), str(waveht[i]),
               horizontalalignment='center', fontsize=13, color='b' )
               #rotation=-90,verticalalignment='center', fontsize=11, color='b' )

        plt.text(xkeys[marks[-1]], ykeys[0]+2.2*(ykeys[1]-ykeys[0]), 'SWH m',
            horizontalalignment='right', fontsize=15, color='k' )
            #rotation=-90,verticalalignment='center', fontsize=15, color='k' )

        # Refresh subplots and save them.
        plt.subplots_adjust(wspace=0.02, hspace=0.01)

        plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', 
                    orientation='landscape',papertype=paprtyp,format='ps')

        plt.close()
    # End of swhglobal plot.


def swhlocal(swhs, verts, ncels, colrs, config,
             mdlname='SMC', datx='2018', psfile='output.ps', 
             paprorn='portrait', paprtyp='a3'):
    """
    ##  Draw local swh plot with given polycollections and cell
    ##  array. Output as A3 ps file.        JGLi28Feb2019
    """

    #  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

    #  Local plot configuration parameters
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]

    #  Maximum mapping radius.
    radius=rdpols[0]

    #  Set up wave height scale and marks with colrs.N.
    #  colrs.N returns colrs' total number of colors 256.
    waveht, factor, residu, marks, ncstr, nclrm = scale_swh(nclrm=colrs.N)

    resmn1=residu - 1.0
    nswh0= ncstr+ int( factor*np.log(-resmn1 + residu) )

    #print ' nswh0 and marks = ', nswh0, marks
    #print ' factor, residu, resmn1 = %f, %f, %f' % (factor, residu, resmn1) 

    #  Some constant variables for plots.
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    
    #  Use ijk to count how many times to draw.
    ijk=0

    if True:
        #  Work out max and min values, excluding missing data (-999.0)
        cmax = swhs.max()
        cmin = swhs[ swhs > -999.0 ].min()
        print ( ' swh range %f, %f' % (cmin, cmax) )
        cmxs = 'SWHmx = %6.2f m' % cmax
        cmns = 'SWHmn = %10.3E' % cmin

        #  Reset missing values (-999.0) to be -resmn1 
        swhs[ swhs < -resmn1] = -resmn1

        #  Trim large values into plot range if any
        swhs[ swhs > 32.0 ] = 32.0 

        #  Convert swhs with logarithm scale.
        icnf = ncstr+ np.rint( factor*np.log(swhs+residu) )
        nswh = np.array( icnf, dtype=np.int16 )

        print (" Drawing "+psfile)

        #  Set up first subplot and axis for northern hemisphere
        fig=plt.figure(figsize=sztpxy[0:2])
        ax=fig.add_subplot(1,1,1)
        ax.set(**xprop)
        ax.set(**yprop)
        ax.set_aspect('equal')
        ax.set_autoscale_on(False)
        ax.set_axis_off()
        plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

        # Prepare PolyCollection for this plot.
        polynorth = PolyCollection(verts)

        # Create color array for this plot
        pcface = []
        pcedge = []
        for i in ncels:
            if( nswh[i] != nswh0 ): 
                pcface.append(colrs(nswh[i]))
                pcedge.append(colrs(nswh[i]))
            else:
                pcface.append(colrs(255)) 
                pcedge.append(colrs(0)) 

        #smcpoly.set_color(pcolr)        # This line defines both edge and face color.
        polynorth.set_facecolor(pcface)
        polynorth.set_edgecolor(pcedge)
        polynorth.set_linewidth( 0.2 )
        #print (" Drawing north hemisphere cells ... ")
        ax.add_collection(polynorth)  

        #  Draw colorbar inside plot.
        xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks)
        ax.add_collection(clrply)
        dkx=clrbxy[2]; dky=clrbxy[3]
        for i in range(len(waveht)):
            m = marks[i]
            if( dkx < dky ): 
                plt.text(xkeys[0]+1.15*dkx, ykeys[m], str(waveht[i]),
                    verticalalignment='center', fontsize=11, color='b' )
            else:
                plt.text(xkeys[m], ykeys[0]+1.15*dky, str(waveht[i]),
                    horizontalalignment='center', fontsize=11, color='b' )

        if( dkx < dky ):
            plt.text(xkeys[0]+2.0*dkx, ykeys[marks[3]], 'SWH m',
                rotation=90,verticalalignment='center', fontsize=15, color='k' )
        else:
            plt.text(xkeys[marks[3]], ykeys[0]+2.0*dky, 'SWH m',
                rotation=0,horizontalalignment='center', fontsize=15, color='k' )

        #  Put cell information inside plot
        tpx= sztpxy[2]
        tpy= sztpxy[3]
        plt.text(tpx, tpy-1.0, mdlname, 
             horizontalalignment='center', fontsize=17, color='g' )
        plt.text(tpx, tpy-1.6, datx,
             horizontalalignment='center', fontsize=15, color='k' )
        plt.text(tpx, tpy-2.1, cmxs,
             horizontalalignment='center', fontsize=13, color='r' )
        plt.text(tpx, tpy-2.6, cmns,
             horizontalalignment='center', fontsize=13, color='b' )

        #  Refresh subplots and save them.
        plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', 
                    orientation=paprorn,papertype=paprtyp)

        plt.close()

    #  End of swhlocal 
