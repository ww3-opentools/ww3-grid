#==================================================================================
# BSD License
#
# Copyright (c)2020, ww3-opentools developers, all rights reserved
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
# gebco_reduce.py
#
# PURPOSE:
#  Functions library for generating a reduced version of the GEBCO bathymetry
#  product, either based on averaging over an integer number of GEBCO grid cells
#  or using scipy interpolate functions. The reduction makes subsequent work for 
#  gridgen quicker and easier, since the cell finding loops in that code then get 
#  to work with smaller source grid files and arrays.
#
#  Functions are also used to correct sea levels and/or remove large inland lakes
#  and areas of land below mean sea level from the source bathymetry, and provide
#  a percentage land value for each of the coarsened cells.
#
# REVISION HISTORY:
#
# A. Saulter; Met Office; May-2020; Version:1.0
#  Code prepared for initial release on github
#
#==================================================================================

import netCDF4 as nc
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

def correctLakesBathy(latout, lonout, depout, mskout, depthmin=0.0, removesmall=None, caspianonly=False):
    '''Apply known height corrections to inland locations'''

    lakecorrections = {}
    # lakes: correct or remove
    lakecorrections['CaspianSea'] = {'lonw':44.0,'lone':55.6,'lats':36.0,'latn':51.0,'correction':27.5,'removal':1000.0}
    if caspianonly:
        lakecorrections['LakeEyre'] = {'lonw':136.0,'lats':-31.9,'lone':140.9,'latn':-26.7,'removal':1000.0}
        lakecorrections['LakeSuperior'] = {'lonw':-92.4,'lats':46.4,'lone':-84.2,'latn':49.1,'removal':1000.0}
        lakecorrections['LakeMichigan'] = {'lonw':-88.15,'lats':41.55,'lone':-84.73,'latn':49.1,'removal':1000.0}
        lakecorrections['LakeHuron'] = {'lonw':-84.73,'lats':42.96,'lone':-79.65,'latn':46.4,'removal':1000.0}
        lakecorrections['LakeErie'] = {'lonw':-83.58,'lats':41.35,'lone':-78.82,'latn':42.94,'removal':1000.0}
        lakecorrections['LakeOntario'] = {'lonw':-79.89,'lats':43.14,'lone':-75.99,'latn':44.34,'removal':1000.0}
        lakecorrections['LakeLadoga'] = {'lonw':30.27,'lats':59.89,'lone':33.02,'latn':51.78,'removal':100.0}
        lakecorrections['LakeOnega'] = {'lonw':33.64,'lats':60.84,'lone':36.54,'latn':62.95,'removal':100.0}
    else:
        lakecorrections['LakeEyre'] = {'lonw':136.0,'lats':-31.9,'lone':140.9,'latn':-26.7,'correction':12,'removal':1000.0}
        lakecorrections['LakeSuperior'] = {'lonw':-92.4,'lats':46.4,'lone':-84.2,'latn':49.1,'correction':-183.0,'removal':1000.0}
        lakecorrections['LakeMichigan'] = {'lonw':-88.15,'lats':41.55,'lone':-84.73,'latn':49.1,'correction':-176.0,'removal':1000.0}
        lakecorrections['LakeHuron'] = {'lonw':-84.73,'lats':42.96,'lone':-79.65,'latn':46.4,'correction':-176.0,'removal':1000.0}
        lakecorrections['LakeErie'] = {'lonw':-83.58,'lats':41.35,'lone':-78.82,'latn':42.94,'correction':-173.0,'removal':1000.0}
        lakecorrections['LakeOntario'] = {'lonw':-79.89,'lats':43.14,'lone':-75.99,'latn':44.34,'correction':-73.0,'removal':1000.0}
        lakecorrections['LakeVictoria'] = {'lonw':31.465,'lats':-3.1,'lone':35.0,'latn':0.5,'correction':-1135.0,'removal':0.0}
        lakecorrections['LakeTanganyika'] = {'lonw':29.00,'lats':-8.88,'lone':31.3,'latn':-3.25,'correction':-773.0,'removal':0.0}
        lakecorrections['LakeBaikal'] = {'lonw':103.56,'lats':51.45,'lone':110.0,'latn':55.95,'correction':-455.0,'removal':0.0}
        lakecorrections['GreatBearLake'] = {'lonw':-125.2,'lats':64.76,'lone':-117.4,'latn':67.07,'correction':-156.0,'removal':1000.0}
        lakecorrections['LakeMalawi'] = {'lonw':33.85,'lats':-14.44,'lone':35.31,'latn':-9.45,'correction':-468.0,'removal':0.0}
        lakecorrections['GreatSlaveLake'] = {'lonw':-117.31,'lats':60.81,'lone':-108.8,'latn':62.94,'correction':-156.0,'removal':1000.0}
        lakecorrections['LakeWinnipeg'] = {'lonw':-99.32,'lats':50.15,'lone':-96.21,'latn':54.09,'correction':-217.0,'removal':1000.0}
        lakecorrections['LakeLadoga'] = {'lonw':30.27,'lats':59.89,'lone':33.02,'latn':51.78,'correction':-5.0,'removal':100.0}
        lakecorrections['LakeOnega'] = {'lonw':33.64,'lats':60.84,'lone':36.54,'latn':62.95,'correction':-33.0,'removal':100.0}
        lakecorrections['LakeBalkash'] = {'lonw':73.32,'lats':44.76,'lone':79.36,'latn':46.92,'correction':-341.4,'removal':0.0}

    # depressions: remove only
    lakecorrections['TurfanDepression'] = {'lonw':87.5,'lats':41.5,'lone':90.7,'latn':43.9,'removal':1000.0} # 154m below sea level
    lakecorrections['CaspianDepression'] = {'lonw':62.3,'lats':41.2,'lone':63.8,'latn':42.8,'removal':1000.0} # lake next to Caspian, 27.5m below sea level
    lakecorrections['ChottMelrhir'] = {'lonw':5.9,'lats':33.8,'lone':7.3,'latn':34.75,'removal':1000.0} # 40m below sea level
    lakecorrections['ShattAlGharsah'] = {'lonw':7.35,'lats':34.0,'lone':8.15,'latn':34.2,'removal':1000.0} # 17m below sea level
    lakecorrections['SabkhatGhuzayyil'] = {'lonw':7.35,'lats':34.0,'lone':8.15,'latn':34.2,'removal':1000.0} # 17m below sea level
    lakecorrections['QattaraDepression'] = {'lonw':23.90,'lats':28.33,'lone':31.58,'latn':30.50,'removal':1000.0} # 133m below sea level
    lakecorrections['DeadSeaDepression'] = {'lonw':35.2,'lats':30.5,'lone':36.1,'latn':32.95,'removal':1000.0} # 413m below sea level
    lakecorrections['SaltonTrough'] = {'lonw':-116.4,'lats':31.85,'lone':-115.1,'latn':33.7,'removal':1000.0} # 69m below sea level
    lakecorrections['DeathValley'] = {'lonw':-117.3,'lats':35.5,'lone':-116.3,'latn':36.65,'removal':1000.0} # 83m below sea level
    lakecorrections['UnknownNorthernCalifornia'] = {'lonw':-122.09,'lats':37.75,'lone':-121.25,'latn':38.50,'removal':1000.0} # unknown
    lakecorrections['StJuliansGreatDepression'] = {'lonw':-68.55,'lats':-49.7,'lone':-68.0,'latn':-49.35,'removal':1000.0} # 105m below sea level
    lakecorrections['GuanBajoDelGualicho'] = {'lonw':-65.75,'lats':-40.6,'lone':-64.7,'latn':-40.00,'removal':1000.0} # 73m below sea level
    lakecorrections['UnknownBuenosAires'] = {'lonw':-64.4,'lats':-38.9,'lone':-62.7,'latn':-38.2,'removal':1000.0} # unknown
    lakecorrections['AmazonDelta'] = {'lonw':-55.0,'lats':-2.5,'lone':-49.5,'latn':-0.35,'removal':1000.0} # unknown
    lakecorrections['AfarDepression'] = {'lonw':42.36,'lats':11.5,'lone':42.68,'latn':11.82,'removal':1000.0} # 155m below sea level
    lakecorrections['DenakilDepression1'] = {'lonw':40.00,'lats':12.65,'lone':40.65,'latn':14.85,'removal':1000.0} # 125m below sea level
    lakecorrections['DenakilDepression2'] = {'lonw':40.65,'lats':12.65,'lone':41.1,'latn':14.36,'removal':1000.0} # 125m below sea level
    lakecorrections['DenakilDepression3'] = {'lonw':41.1,'lats':12.65,'lone':41.5,'latn':13.7,'removal':1000.0} # 125m below sea level
    lakecorrections['UnknownKuwait'] = {'lonw':46.6,'lats':30.34,'lone':48.45,'latn':31.1,'removal':1000.0} # unknown
    lakecorrections['UnknownSaudiArabia1'] = {'lonw':51.3,'lats':23.7,'lone':51.6,'latn':23.95,'removal':1000.0} # unknown
    lakecorrections['UnknownSaudiArabia2'] = {'lonw':50.85,'lats':24.45,'lone':51.05,'latn':24.65,'removal':1000.0} # unknown
    lakecorrections['DonRiver'] = {'lonw':39.7,'lats':47.15,'lone':40.35,'latn':47.5,'removal':1000.0} # unknown
    lakecorrections['UnknownKazakhstan'] = {'lonw':54.76,'lats':45.4,'lone':56.57,'latn':46.5,'removal':1000.0} # unknown
    lakecorrections['UnknownTurkmenistan'] = {'lonw':54.2,'lats':42.5,'lone':54.9,'latn':43.3,'removal':1000.0} # unknown
    lakecorrections['AkdzhakayaDepression'] = {'lonw':55.2,'lats':40.0,'lone':58.3,'latn':41.15,'removal':1000.0} # 81m below sea level; wider area 
    lakecorrections['LakeTaiYangtze'] = {'lonw':119.0,'lats':30.92,'lone':120.75,'latn':32.9,'removal':1000.0} # unknown
    lakecorrections['Yancheng1'] = {'lonw':119.1,'lats':32.9,'lone':120.25,'latn':34.27,'removal':1000.0} # unknown
    lakecorrections['Yancheng2'] = {'lonw':120.25,'lats':32.9,'lone':120.5,'latn':33.9,'removal':1000.0} # unknown
    lakecorrections['KaragiyeDepression1'] = {'lonw':51.38,'lats':43.37,'lone':51.98,'latn':43.9,'removal':1000.0} # 132m below sea level
    lakecorrections['KaragiyeDepression2'] = {'lonw':51.67,'lats':43.21,'lone':51.98,'latn':43.37,'removal':1000.0} # 132m below sea level
    lakecorrections['SzczecinLagoon'] = {'lonw':13.85,'lats':53.6,'lone':14.63,'latn':53.92,'removal':1000.0} # unknown
    lakecorrections['Halligen'] = {'lonw':8.58,'lats':54.25,'lone':9.0,'latn':55.10,'removal':1000.0} # unknown
    lakecorrections['Elbe'] = {'lonw':9.06,'lats':53.7,'lone':9.9,'latn':53.9,'removal':1000.0} # unknown
    lakecorrections['UnknownNLGer1'] = {'lonw':7.6,'lats':53.4,'lone':10.2,'latn':53.7,'removal':1000.0} # unknown
    lakecorrections['UnknownNLGer2'] = {'lonw':5.53,'lats':59.2,'lone':9.3,'latn':53.4,'removal':1000.0} # unknown
    lakecorrections['EastAnglia'] = {'lonw':-0.5,'lats':52.34,'lone':0.6,'latn':52.7,'removal':1000.0} # unknown

    for region in lakecorrections.keys():
        correctlake = True
        if lakecorrections[region]['lone']<np.min(lonout) or lakecorrections[region]['latn']<np.min(latout) or \
           lakecorrections[region]['lonw']>np.max(lonout) or lakecorrections[region]['lats']>np.max(latout):
           correctlake=False
        if correctlake:
            print('[INFO] Correcting heights/depths for %s' %region)
            llx = np.min(np.where(lonout>lakecorrections[region]['lonw']))
            urx = np.max(np.where(lonout<lakecorrections[region]['lone']))+1
            lly = np.min(np.where(latout>lakecorrections[region]['lats']))
            ury = np.max(np.where(latout<lakecorrections[region]['latn']))+1
            if 'correction' in lakecorrections[region].keys():
                depout[lly:ury,llx:urx] = depout[lly:ury,llx:urx] + lakecorrections[region]['correction']
                if lakecorrections[region]['correction'] < 0.0:
                    mskout[depout <= depthmin] = 0.0
                else:
                    mskout[depout > depthmin] = 1.0
            else:
                depout[lly:ury,llx:urx] = depout[lly:ury,llx:urx] + lakecorrections[region]['removal']
                mskout[depout > depthmin] = 1.0

    if removesmall is not None:
        print('[INFO] Checking through grid for small isolated water bodies at grid size %d' %removesmall)
        chksize = removesmall + 2
        delcounter = 0
        for lpy in range(np.shape(mskout)[0] - chksize): 
            if np.mod(lpy,100) == 0:
                print('[INFO] ..processed rows for %d y-cells and removed %d small water bodies' %(lpy,delcounter))
            for lpx in range(np.shape(mskout)[1] - chksize): 
                chktmp = np.copy(mskout[lpy:lpy+chksize,lpx:lpx+chksize])
                if (not np.all(chktmp == 1.0)) and (not np.all(chktmp == 0.0)): 
                    chktmp[1:1+removesmall,1:1+removesmall] = 1.0
                    if np.all(chktmp == 1.0):
                        mskout[lpy+1:lpy+1+removesmall,lpx+1:lpx+1+removesmall] = 1.0
                        #depout[lpy+1:lpy+1+removesmall,lpx+1:lpx+1+removesmall] = 10.0
                        delcounter = delcounter + 1
        print('[INFO] Removed %d small water bodies from land-sea mask' %delcounter)

    return depout, mskout


def correctLakesBathyfromfile(ncfile, datadir='.', depthmin=0.0, removesmall=2, caspianonly=True, pltchk=True):
    '''Apply inland height corrections direct to file'''

    print('[INFO] Applying inland height/depth corrections direct to netCDF file')
    print('[WARN] This action will overwrite existing depth and land-sea mask data')
    d = nc.Dataset(datadir + '/' + ncfile,'a')
    lonout = d.variables['lon'][:]
    latout = d.variables['lat'][:]
    depout = d.variables['elevation'][:,:]
    mskout = d.variables['landmask'][:,:]
    depout, mskout = correctLakesBathy(latout, lonout, depout, mskout, 
                                        depthmin=depthmin, removesmall=removesmall, caspianonly=caspianonly)
    elev = d.variables['elevation']
    mask = d.variables['landmask']
    elev[:,:] = depout[:,:]
    mask[:,:] = mskout[:,:]
    if pltchk:
        plotGrid(latout, lonout, depths=depout, landsea=mskout, depthmin=5.0, depthmax=-500.0)

    d.corrected_land_lakes = 'True'
    d.close()

    return


def writeReducedNC(outfile, scalefac, depthmin, latout, lonout, 
                    depout, mskout, datadir='.'):
    '''Write out the reduced dataset to a new netCDF file'''

    loresfile = datadir + '/' + outfile
    print('[INFO] Writing reduced grid data to %s' %loresfile)
    with nc.Dataset(loresfile, 'w') as nbg:
        ndimx = nbg.createDimension('lon',size=np.size(lonout))
        ndimy = nbg.createDimension('lat',size=np.size(latout))

        ndep = nbg.createVariable('lat','f8',dimensions=('lat'))
        ndep.units = 'degrees_east'
        ndep[:] = latout[:]

        ndep = nbg.createVariable('lon','f8',dimensions=('lon'))
        ndep.units = 'degrees_north'
        ndep[:] = lonout[:]

        ndep = nbg.createVariable('elevation','f8',dimensions=('lat','lon'))
        ndep.units = 'm'
        ndep[:,:] = depout[:,:]

        ndep = nbg.createVariable('landmask','f8',dimensions=('lat','lon'))
        ndep.units = '1'
        ndep[:,:] = mskout[:,:]

        # add global attributes to describe processing
        nbg.description = 'Reduced GEBCO bathymetry grid: mean depth values over cell'
        nbg.reduction_scale_factor = scalefac
        nbg.minimum_depth = depthmin

        #nbg.close()

    return


def writeInterpolatedNC(outfile, dx, dy, depthmin, latout, lonout, 
                    depout, mskout, datadir='.'):
    '''Write out the reduced dataset to a new netCDF file'''

    loresfile = datadir + '/' + outfile
    print('[INFO] Writing reduced grid data to %s' %loresfile)
    with nc.Dataset(loresfile, 'w') as nbg:
        ndimx = nbg.createDimension('lon',size=np.size(lonout))
        ndimy = nbg.createDimension('lat',size=np.size(latout))

        ndep = nbg.createVariable('lat','f8',dimensions=('lat'))
        ndep.units = 'degrees_east'
        ndep[:] = latout[:]

        ndep = nbg.createVariable('lon','f8',dimensions=('lon'))
        ndep.units = 'degrees_north'
        ndep[:] = lonout[:]

        ndep = nbg.createVariable('elevation','f8',dimensions=('lat','lon'))
        ndep.units = 'm'
        ndep[:,:] = depout[:,:]

        ndep = nbg.createVariable('landmask','f8',dimensions=('lat','lon'))
        ndep.units = '1'
        ndep[:,:] = mskout[:,:]

        # add global attributes to describe processing
        nbg.description = 'Interpolated GEBCO bathymetry grid: mean depth values over cell'
        nbg.interpolation_dx = dx
        nbg.interpolation_dy = dy
        nbg.minimum_depth = depthmin

        #nbg.close()

    return


def plotGrid(latout, lonout, depths=None, landsea=None, depthmin=5.0, depthmax=-500.0):
    '''Plots out the reduced grid as a check'''

    print('[INFO] Plotting depths and/or land-sea mask for checking..')
    if (depths is not None) and (landsea is not None):
        plt.subplot(1,2,1)
    if depths is not None:
        depouttmp = np.ma.masked_greater(depths, depthmin)
        plt.pcolormesh(lonout, latout, depouttmp, vmin=depthmax)
        plt.colorbar()
    if (depths is not None) and (landsea is not None):
        plt.subplot(1,2,2)
    if landsea is not None:
        plt.pcolormesh(lonout, latout, landsea)
        plt.colorbar()
    plt.show()

    return


def plotGridfromfile(ncfile, datadir='.', usedepths=True, uselandsea=True, depthmin=5.0, depthmax=-500.0):
    '''Plots out the reduced grid as a check'''

    depout = None
    mskout = None
    myncfile = datadir + '/' + ncfile
    print('[INFO] Getting depths and/or land-sea mask data from %s' %myncfile)
    d = nc.Dataset(myncfile)
    lonout = d.variables['lon'][:]
    latout = d.variables['lat'][:]
    if usedepths:
        depout = d.variables['elevation'][:,:]
    if uselandsea:
        mskout = d.variables['landmask'][:,:]
    plotGrid(latout, lonout, depths=depout, landsea=mskout, depthmin=depthmin, depthmax=depthmax)
    d.close()

    return


def rebin(arr, new_shape):
    '''Bin a large array to a smaller one by averaging'''

    if np.size(new_shape) == 1:
        shape = (new_shape[0], arr.shape[0] // new_shape[0])
        arrout = arr.reshape(shape).mean(-1)
    elif np.size(new_shape) == 2:
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
                 new_shape[1], arr.shape[1] // new_shape[1])
        arrout = arr.reshape(shape).mean(-1).mean(1)

    return arrout


def reduceGrid(hiresfile, scalefac=6, depthmin=0.0, cutout=None, loopmethod=False):
    '''Reduce the size of GEBCO grid by averaging cells according to integer scale factor.
       Presently this has only been tested with GEBCO_2014, i.e. 30 seconds grid'''

    print('[INFO] Running reduction of GEBCO data at scale factor %d' %scalefac)
    print('[INFO] Reading data from %s' %hiresfile)
    d = nc.Dataset(hiresfile)
    nlats = d.dimensions['lat'].size
    nlons = d.dimensions['lon'].size
    dlat = 180.0 / np.float(nlats)
    dlon = 360.0 / np.float(nlons)

    tilemaxx = 10800
    tilemaxy = 5400
    if cutout is None:
        dimy = nlats
        dimx = nlons
        offsx = 0
        offsy = 0
    else:
        if cutout[0] > 180.0: cutout[0] = cutout[0]-360.0
        if cutout[2] > 180.0: cutout[0] = cutout[0]-360.0
        offsx = np.int(np.floor((cutout[0]+180.0)/dlon))-1
        offsy = np.int(np.floor((cutout[1]+90.0)/dlat))-1
        dimx = np.int(np.ceil((cutout[2]-cutout[0])/dlon))
        dimx = np.int((np.ceil(dimx/scalefac)+1)*scalefac)
        dimy = np.int(np.ceil((cutout[3]-cutout[1])/dlat))
        dimy = np.int((np.ceil(dimy/scalefac)+1)*scalefac)

    # set the output arrays for reduced grid
    newy = np.int(dimy/scalefac)
    newx = np.int(dimx/scalefac)
    depout = np.zeros([newy, newx])
    mskout = np.zeros([newy, newx])
    latout = np.zeros(newy)
    lonout = np.zeros(newx)

    # using loops/tiles here to avoid memory problems with read in of full dataset
    if loopmethod:
        # loop using strides as the averaging routine
        lptot = np.float(scalefac**2)
        for lpx in range(scalefac):
            for lpy in range(scalefac):
                lpiter = lpx * scalefac + lpy
                print('[INFO] Working on iteration %d' %lpiter + ' out of %d' %lptot)
                dep = d.variables['elevation'][lpy+offsy:dimy+offsy:scalefac,lpx+offsx:dimx+offsx:scalefac]
                msk = np.zeros(np.shape(dep))
                msk[dep > depthmin] = 1.
                lat = d.variables['lat'][lpy+offsy:dimy+offsy:scalefac]
                lon = d.variables['lon'][lpx+offsx:dimx+offsx:scalefac]
                depout = depout + dep / lptot
                mskout = mskout + msk / lptot
                lonout = lonout + lon / lptot
                latout = latout + lat / lptot
    else:
        # run the tiling method
        lpx = offsx
        lpoutx = 0
        while lpx < offsx+dimx:
            addx = np.min([tilemaxx, offsx+dimx-lpx])
            addoutx = np.int(addx/scalefac)
            lpy = offsy
            lpouty = 0        
            while lpy < offsy+dimy:
                print('[INFO] Working on tile from x:%d, y:%d' %(lpx,lpy))
                addy = np.min([tilemaxy, offsy+dimy-lpy])
                addouty = np.int(addy/scalefac)
                tmpdep = d.variables['elevation'][lpy:lpy+addy,lpx:lpx+addx]
                tmpmsk = np.zeros(np.shape(tmpdep))
                tmpmsk[tmpdep > depthmin] = 1.
                tmplat = d.variables['lat'][lpy:lpy+addy]
                tmplon = d.variables['lon'][lpx:lpx+addx]
                print('... read subdomain from file')
                depout[lpouty:lpouty+addouty,lpoutx:lpoutx+addoutx] = rebin(tmpdep, [addouty,addoutx])
                mskout[lpouty:lpouty+addouty,lpoutx:lpoutx+addoutx] = rebin(tmpmsk, [addouty,addoutx])
                latout[lpouty:lpouty+addouty] = rebin(tmplat, [addouty])
                lonout[lpoutx:lpoutx+addoutx] = rebin(tmplon, [addoutx])
                print('... reduced subdomain')
                lpy = lpy + addy
                lpouty = lpouty + addouty
            lpx = lpx + addx
            lpoutx = lpoutx + addoutx

    d.close()

    return latout, lonout, depout, mskout


def reduceGEBCO(scalefac=6, depthmin=0.0, cutout=None, region=None,
                 pltchk=True, correctlakes=False,
                 gebcofile='GEBCO_2014_2D.nc', datadir='.'):
    '''Controls the grid reduction process'''

    # generate the reduced grid
    hiresfile = datadir + '/' + gebcofile
    latout, lonout, depout, mskout = reduceGrid(hiresfile, scalefac=scalefac, depthmin=depthmin, cutout=cutout)
    if pltchk:
        plotGrid(latout, lonout, depths=depout, landsea=mskout, depthmin=5.0, depthmax=-500.0)

    # correct the land bathy as required
    # at present this will run with the correctLakesBathy defaults 
    #  - adds large lakes and does not remove small water bodies
    if correctlakes:
        depout, mskout = correctLakesBathy(latout, lonout, depout, mskout)
        if pltchk:
            plotGrid(latout, lonout, depths=depout, landsea=mskout, depthmin=5.0, depthmax=-500.0)

    # write out data to a new netCDF file
    if region is None:
        loresfile = 'GEBCO_reduced_%d' %scalefac + '.nc'
    else:
        loresfile = 'GEBCO_reduced_%s' %region + '_%d' %scalefac + '.nc'
    writeReducedNC(loresfile, scalefac, depthmin, latout, lonout, 
                    depout, mskout, datadir=datadir)    

    return


def interpGrid(hiresfile, dx=0.5, dy=0.5, depthmin=0.0, cutout=None, loopmethod=False):
    '''Reduce the size of GEBCO grid by interpolating cells'''

    print('[INFO] Running reduction of GEBCO data to resolution dx:%.4f, dy:%.4f' %(dx,dy))
    print('[INFO] Reading data from %s' %hiresfile)
    d = nc.Dataset(hiresfile)
    nlats = d.dimensions['lat'].size
    nlons = d.dimensions['lon'].size
    dlat = 180.0 / np.float(nlats)
    dlon = 360.0 / np.float(nlons)

    tilemaxx = 7200
    tilemaxy = 3600
    if cutout is None:
        x0 = -180.0 + dx / 2.0
        y0 = -90.0 + dy / 2.0
        xl = -1.0 * x0
        yl = -1.0 * y0
        offsx = 0
        offsy = 0
    else:
        if cutout[0] > 180.0: cutout[0] = cutout[0]-360.0
        if cutout[2] > 180.0: cutout[0] = cutout[0]-360.0
        x0 = cutout[0]
        y0 = cutout[1]
        xl = cutout[2]
        yl = cutout[3]
        offsx = np.int(np.floor(cutout[0]/dlon))-1
        offsy = np.int(np.floor(cutout[1]/dlat))-1

    # set the output arrays for reduced grid
    newy = np.int((yl-y0)/dy+1)
    newx = np.int((xl-x0)/dx+1)
    depout = np.zeros([newy, newx])
    mskout = np.zeros([newy, newx])
    latout = np.arange(y0, yl+dy/5.0, dy)
    lonout = np.arange(x0, xl+dx/5.0, dx)

    # using loops/tiles here to avoid memory problems with read in of full dataset
    lpx = offsx
    lpoutx = 0
    while lpx < nlons-1:
        addx = np.min([tilemaxx+1, nlons-lpx])
        tmplon = d.variables['lon'][lpx:lpx+addx]
        addoutx = np.int(np.floor((tmplon[-1] - lonout[lpoutx]) / dx))
        print('[INFO] Working on tile from x:%d' %lpx)
        print('[INFO] GEBCO longitude range: %.8f to %.8f' %(tmplon[0],tmplon[-1]))
        print('[INFO] Number of processed grid cells in x dimension: %d' %addoutx)
        print('[INFO] Processed grid longitude range: %.8f to %.8f' 
              %(lonout[lpoutx],lonout[lpoutx+addoutx]))
        lpy = offsy
        lpouty = 0        
        while lpy < nlats-1:
            print('[INFO] Working on tile from x:%d, y:%d' %(lpx,lpy))
            addy = np.min([tilemaxy+1, nlons-lpy])
            tmplat = d.variables['lat'][lpy:lpy+addy]
            addouty = np.int(np.floor((tmplat[-1] - latout[lpouty]) / dy))
            print('[INFO] GEBCO latitude range: %.8f to %.8f' %(tmplat[0],tmplat[-1]))
            print('[INFO] Number of processed grid cells in y dimension: %d' %addouty)
            print('[INFO] Processed grid longitude range: %.8f to %.8f' 
                  %(latout[lpouty],latout[lpouty+addouty]))
            tmpdep = d.variables['elevation'][lpy:lpy+addy,lpx:lpx+addx]
            tmpmsk = np.zeros(np.shape(tmpdep))
            tmpmsk[tmpdep > depthmin] = 1.
            print('[INFO]... read subdomain from file')
            splinedep = interp.RectBivariateSpline(tmplat,tmplon,tmpdep)
            depout[lpouty:lpouty+addouty+1,lpoutx:lpoutx+addoutx+1] = \
              splinedep(latout[lpouty:lpouty+addouty+1],lonout[lpoutx:lpoutx+addoutx+1])
            splinemsk = interp.RectBivariateSpline(tmplat,tmplon,tmpmsk)
            mskout[lpouty:lpouty+addouty+1,lpoutx:lpoutx+addoutx+1] = \
              splinemsk(latout[lpouty:lpouty+addouty+1],lonout[lpoutx:lpoutx+addoutx+1])
            print('[INFO]... interpolated to reduced subdomain')
            lpy = lpy + addy - 1
            lpouty = lpouty + addouty + 1
        lpx = lpx + addx - 1
        lpoutx = lpoutx + addoutx + 1

    # correct spline interpolation limits for mask
    mskout[mskout < 0.005] = 0.0
    mskout[mskout > 0.995] = 1.0

    d.close()

    return latout, lonout, depout, mskout


def interpGEBCO(dx=0.5, dy=0.5, depthmin=0.0, region=None, cutout=None,
                 pltchk=True, correctlakes=False,
                 gebcofile='GEBCO_2014_2D.nc', datadir='.'):
    '''Controls the grid interpolation process'''

    # generate the reduced grid
    hiresfile = datadir + '/' + gebcofile
    latout, lonout, depout, mskout = interpGrid(hiresfile, dx=dx, dy=dy, depthmin=depthmin, cutout=cutout)
    if pltchk:
        plotGrid(latout, lonout, depths=depout, landsea=mskout, depthmin=5.0, depthmax=-500.0)

    # correct the land bathy as required
    # at present this will run with the correctLakesBathy defaults 
    #  - adds large lakes and does not remove small water bodies
    if correctlakes:
        depout, mskout = correctLakesBathy(latout, lonout, depout, mskout)
        if pltchk:
            plotGrid(latout, lonout, depths=depout, landsea=mskout, depthmin=5.0, depthmax=-500.0)

    # write out data to a new netCDF file
    interpstr = '%d' %np.floor(dx) + 'd%04d' %np.floor(10000*(dx-np.floor(dx)))
    interpstr = interpstr + 'by%d' %np.floor(dy) + 'd%04d' %np.floor(10000*(dy-np.floor(dy)))
    if region is not None:
        interpstr = region + '_' + interpstr
    loresfile = 'GEBCO_interpolated_%s' %interpstr + '.nc'
    writeInterpolatedNC(loresfile, dx, dy, depthmin, latout, lonout, 
                    depout, mskout, datadir=datadir)    

    return
