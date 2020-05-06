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
# run_gridgen.py
#
# PURPOSE:
#  Runs gridgen.py functions based on user defined configuration file inputs
#
# REVISION HISTORY:
#
# A. Saulter; Met Office; May-2020; Version: 1.0
#  Code prepared for initial release on github
#
#==================================================================================

import matplotlib.pyplot as plt
import numpy as np
import configparser as cfg
import sys

import gridgen as grd


def readConfigMark(config, procname):
    '''Reads in relevant content of config file for cell marking actions
       This is called as a separate function so that mark actions can be
       composited'''

    # place mark actions into a list
    cfginfo = {'action':'mark',
               'marklist':[]}
    if config.get(procname,'action').lower() == 'markmulti':
        for markproc in config.get(procname,'marknames').split(','):
            cfginfo['marklist'].append(markproc.replace(' ',''))
    else:
        cfginfo['marklist'].append(procname)

    # run through list and append actions according to type
    for proc in cfginfo['marklist']:

        # mark cells below a given depth
        if config.get(proc,'action').lower() == 'markdepths':
            # set defaults
            mrkinfo = {'action':'markdepths',
                       'workdir':'.',
                       'basefile':'smcGrid_basegrid.nc',
                       'name':'smcGrid',
                       'label':'newmark',
                       'markertype':'dry',
                       'depthlim':5.0}
            # set user options
            if config.has_option(proc,'workdir'):
                mrkinfo['workdir'] = config.get(proc,'workdir')
            if config.has_option(proc,'basefile'):
                mrkinfo['basefile'] = config.get(proc,'basefile')
            if config.has_option(proc,'name'):
                mrkinfo['name'] = config.get(proc,'name')
            if config.has_option(proc,'label'):
                mrkinfo['label'] = config.get(proc,'label')
            if config.has_option(proc,'markertype'):
                mrkinfo['markertype'] = config.get(proc,'markertype')
            if config.has_option(proc,'depthlim'):
                mrkinfo['depthlim'] = np.float(config.get(proc,'depthlim'))

        # mark cells for a specific region (and depth)
        if config.get(proc,'action').lower() == 'markregion':
            # set defaults
            mrkinfo = {'action':'markregion',
                       'workdir':'.',
                       'basefile':'smcGrid_basegrid.nc',
                       'name':'smcGrid',
                       'label':'newmark',
                       'markertype':'tier',
                       'extents':[0.0,-60.0,360.0,73.0],
                       'depthlim':None}
            # set user options
            if config.has_option(proc,'workdir'):
                mrkinfo['workdir'] = config.get(proc,'workdir')
            if config.has_option(proc,'basefile'):
                mrkinfo['basefile'] = config.get(proc,'basefile')
            if config.has_option(proc,'name'):
                mrkinfo['name'] = config.get(proc,'name')
            if config.has_option(proc,'label'):
                mrkinfo['label'] = config.get(proc,'label')
            if config.has_option(proc,'markertype'):
                mrkinfo['markertype'] = config.get(proc,'markertype')
            if config.has_option(proc,'extents'):
                extentstr = config.get(proc,'extents').split(',')
                mrkinfo['extents'] = [np.float(extentstr[0]), np.float(extentstr[1]),
                                      np.float(extentstr[2]), np.float(extentstr[3])]
            if config.has_option(proc,'depthlim'):
                if config.get(proc,'depthlim').lower() != 'none':
                    mrkinfo['depthlim'] = np.float(config.get(proc,'depthlim'))

        # unmark tier cells (set to wet/dry) for a specific region
        if config.get(proc,'action').lower() == 'unmark':
            # set defaults
            mrkinfo = {'action':'unmark',
                       'workdir':'.',
                       'basefile':'smcGrid_basegrid.nc',
                       'name':'smcGrid',
                       'label':'newmark',
                       'markertype':'tier',
                       'extents':[0.0,-60.0,360.0,73.0],
                       'osbox':False,
                       'thruzero':False,
                       'deldry':False}
            # set user options
            if config.has_option(proc,'workdir'):
                mrkinfo['workdir'] = config.get(proc,'workdir')
            if config.has_option(proc,'basefile'):
                mrkinfo['basefile'] = config.get(proc,'basefile')
            if config.has_option(proc,'name'):
                mrkinfo['name'] = config.get(proc,'name')
            if config.has_option(proc,'label'):
                mrkinfo['label'] = config.get(proc,'label')
            if config.has_option(proc,'markertype'):
                mrkinfo['markertype'] = config.get(proc,'markertype')
            if config.has_option(proc,'extents'):
                extentstr = config.get(proc,'extents').split(',')
                mrkinfo['extents'] = [np.float(extentstr[0]), np.float(extentstr[1]),
                                      np.float(extentstr[2]), np.float(extentstr[3])]
            if config.has_option(proc,'osbox'):
                print(config.get(proc,'osbox'))
                if config.get(proc,'osbox').lower() == 'true':
                    mrkinfo['osbox'] = True
            if config.has_option(proc,'thruzero'):
                if config.get(proc,'thruzero').lower() == 'true':
                    mrkinfo['thruzero'] = True
            if config.has_option(proc,'deldry'):
                if config.get(proc,'deldry').lower() == 'true':
                    mrkinfo['deldry'] = True

        cfginfo[proc] = mrkinfo

    return cfginfo


def readConfig(procname, cfgfile):
    '''Reads in relevant content of config file according to action'''

    # read info from the configuration file
    print('[INFO] Searching for configuration information from: %s' %cfgfile)
    config = cfg.RawConfigParser()
    config.read(cfgfile)

    # set SMC base grid
    if config.get(procname,'action').lower() == 'smcbase':
        # set defaults
        cfginfo = {'action':'smcbase',
                   'workdir':'.',
                   'bathyfile':'gebco_reduced.nc',
                   'extents':[0.0,-90.0,360.0,90.0],
                   'dx':0.25,
                   'dy':0.25, 
                   'name':'smcGrid',
                   'label':'basegrid',
                   'mindepth':None, 
                   'drydepthlim':None,
                   'drypcmin':None,
                   'drypcmax':None, 
                   'bathytype':'gebco',
                   'getpcland':True, 
                   'setadj':True}
        # set user options
        if config.has_option(procname,'workdir'):
            cfginfo['workdir'] = config.get(procname,'workdir')
        if config.has_option(procname,'bathyfile'):
            cfginfo['bathyfile'] = config.get(procname,'bathyfile')
        if config.has_option(procname,'extents'):
            extentstr = config.get(procname,'extents').split(',')
            cfginfo['extents'] = [np.float(extentstr[0]), np.float(extentstr[1]),
                                  np.float(extentstr[2]), np.float(extentstr[3])]
        if config.has_option(procname,'dx'):
            cfginfo['dx'] = np.float(config.get(procname,'dx'))
        if config.has_option(procname,'dy'):
            cfginfo['dy'] = np.float(config.get(procname,'dy'))
        if config.has_option(procname,'mindepth'):
            if config.get(procname,'mindepth').lower() != 'none':
                cfginfo['mindepth'] = np.float(config.get(procname,'mindepth'))
        if config.has_option(procname,'drydepthlim'):
            if config.get(procname,'drydepthlim').lower() != 'none':
                cfginfo['drydepthlim'] = np.float(config.get(procname,'drydepthlim'))
        if config.has_option(procname,'drypcmin'):
            if config.get(procname,'drypcmin').lower() != 'none':
                cfginfo['drypcmin'] = np.float(config.get(procname,'drypcmin'))
        if config.has_option(procname,'drypcmax'):
            if config.get(procname,'drypcmax').lower() != 'none':
                cfginfo['drypcmax'] = np.float(config.get(procname,'drypcmax'))
        if config.has_option(procname,'name'):
            cfginfo['name'] = config.get(procname,'name')
        if config.has_option(procname,'label'):
            cfginfo['label'] = config.get(procname,'label')
        if config.has_option(procname,'bathytype'):
            cfginfo['bathytype'] = config.get(procname,'bathytype')
        if config.has_option(procname,'getpcland'):
            if config.get(procname,'getpcland').lower() == 'false':
                cfginfo['getpcland'] = False
        if config.has_option(procname,'setadj'):
            if config.get(procname,'setadj').lower() == 'false':
                cfginfo['setadj'] = False

    # generate an SMC tier
    if config.get(procname,'action').lower() == 'tiergen':
        # set defaults
        cfginfo = {'action':'tiergen',
                   'workdir':'.',
                   'basefile':'smcGrid_basegrid.nc',
                   'bathyfile':'gebco_reduced.nc',
                   'name':'smcGrid',
                   'label':'newtier',
                   'mindepth':None, 
                   'drydepthlim':None,
                   'drypcmin':None,
                   'drypcmax':None, 
                   'bathytype':'gebco',
                   'getpcland':True, 
                   'setadj':True,
                   'deldry':False}
        # set user options
        if config.has_option(procname,'workdir'):
            cfginfo['workdir'] = config.get(procname,'workdir')
        if config.has_option(procname,'basefile'):
            cfginfo['basefile'] = config.get(procname,'basefile')
        if config.has_option(procname,'bathyfile'):
            cfginfo['bathyfile'] = config.get(procname,'bathyfile')
        if config.has_option(procname,'mindepth'):
            if config.get(procname,'mindepth').lower() != 'none':
                cfginfo['mindepth'] = np.float(config.get(procname,'mindepth'))
        if config.has_option(procname,'drydepthlim'):
            if config.get(procname,'drydepthlim').lower() != 'none':
                cfginfo['drydepthlim'] = np.float(config.get(procname,'drydepthlim'))
        if config.has_option(procname,'drypcmin'):
            if config.get(procname,'drypcmin').lower() != 'none':
                cfginfo['drypcmin'] = np.float(config.get(procname,'drypcmin'))
        if config.has_option(procname,'drypcmax'):
            if config.get(procname,'drypcmax').lower() != 'none':
                cfginfo['drypcmax'] = np.float(config.get(procname,'drypcmax'))
        if config.has_option(procname,'name'):
            cfginfo['name'] = config.get(procname,'name')
        if config.has_option(procname,'label'):
            cfginfo['label'] = config.get(procname,'label')
        if config.has_option(procname,'bathytype'):
            cfginfo['bathytype'] = config.get(procname,'bathytype')
        if config.has_option(procname,'getpcland'):
            if config.get(procname,'getpcland').lower() == 'false':
                cfginfo['getpcland'] = False
        if config.has_option(procname,'setadj'):
            if config.get(procname,'setadj').lower() == 'false':
                cfginfo['setadj'] = False
        if config.has_option(procname,'deldry'):
            if config.get(procname,'deldry').lower() == 'true':
                cfginfo['setadj'] = True

    # combine SMC base and tier files
    if config.get(procname,'action').lower() == 'tiercombine':
        # set defaults
        cfginfo = {'action':'tiercombine',
                   'workdir':'.',
                   'basefile':'smcGrid_basegrid.nc',
                   'tierfile':'smcGrid_newtier.nc',
                   'bathyfile':'gebco_reduced.nc',
                   'name':'smcGrid',
                   'label':'newtier',
                   'tiernext':True}
        # set user options
        if config.has_option(procname,'workdir'):
            cfginfo['workdir'] = config.get(procname,'workdir')
        if config.has_option(procname,'basefile'):
            cfginfo['basefile'] = config.get(procname,'basefile')
        if config.has_option(procname,'tierfile'):
            cfginfo['tierfile'] = config.get(procname,'tierfile')
        if config.has_option(procname,'bathyfile'):
            cfginfo['bathyfile'] = config.get(procname,'bathyfile')
        if config.has_option(procname,'name'):
            cfginfo['name'] = config.get(procname,'name')
        if config.has_option(procname,'label'):
            cfginfo['label'] = config.get(procname,'label')
        if config.has_option(procname,'tiernext'):
            if config.get(procname,'tiernext').lower() == 'false':
                cfginfo['tiernext'] = False

    # write the grid in WW3 format
    if config.get(procname,'action').lower() == 'writeww3':
        # set defaults
        cfginfo = {'action':'writeWW3',
                   'workdir':'.',
                   'gridfile':'smcGrid_newgrid.nc',
                   'writedir':'.',
                   'mindepth':None,
                   'writemindepth':False,
                   'arctic':False,
                   'arclat':86.4}
        # set user options
        if config.has_option(procname,'workdir'):
            cfginfo['workdir'] = config.get(procname,'workdir')
        if config.has_option(procname,'gridfile'):
            cfginfo['gridfile'] = config.get(procname,'gridfile')
        if config.has_option(procname,'writedir'):
            cfginfo['writedir'] = config.get(procname,'writedir')
        if config.has_option(procname,'mindepth'):
            if config.get(procname,'mindepth').lower() != 'none':
                cfginfo['mindepth'] = np.float(config.get(procname,'mindepth'))
        if config.has_option(procname,'writemindepth'):
            if config.get(procname,'writemindepth').lower() == 'true':
                cfginfo['writemindepth'] = True
        if config.has_option(procname,'arctic'):
            if config.get(procname,'arctic').lower() == 'true':
                cfginfo['arctic'] = True
        if config.has_option(procname,'arclat'):
            cfginfo['arclat'] = np.float(config.get(procname,'arclat'))

    # cell mark actions
    marklist = ['markmulti','markdepths','markregion','marklandpc','unmark']
    if config.get(procname,'action').lower() in marklist:
        cfginfo = readConfigMark(config, procname)

    return cfginfo        


def runGridgen(cfginfo):
    '''Select and run gridgen actions based on .cfg file namelist'''

    if cfginfo['action'] == 'smcbase':
        # generate the base grid
        print('')
        print('*** Creating SMC base grid')
        basesmc = grd.createBasesmc(cfginfo['bathyfile'],
                                    cfginfo['extents'], 
                                    cfginfo['dx'], 
                                    cfginfo['dy'], 
                                    name=cfginfo['name'],
                                    label=cfginfo['label'],
                                    mindepth=cfginfo['mindepth'],
                                    dlim=cfginfo['drydepthlim'], 
                                    drymin=cfginfo['drypcmin'], 
                                    drymax=cfginfo['drypcmax'], 
                                    bathytype=cfginfo['bathytype'],
                                    getpland=cfginfo['getpcland'],
                                    setadj=cfginfo['setadj'])
        # write base grid file
        basefile = basesmc.writeNC(writedir=cfginfo['workdir'])
        # visualize the grid
        grd.plotGridsmc(basesmc, latlon=False)

    elif cfginfo['action'] == 'mark':
        # mark actions run through a loop to allow composite actions
        for lp, proc in enumerate(cfginfo['marklist']):
            if lp == 0:
                # load the parent grid
                print('')
                print('*** Loading SMC base grid')
                basefile = cfginfo[proc]['workdir'] + '/' + cfginfo[proc]['basefile']
                basesmc = grd.loadNCsmc(basefile)
                print('')
                print('*** Marking cells for SMC grid')
            # marking actions
            if cfginfo[proc]['action'] == 'markdepths':
                basesmc.markDepths(cfginfo[proc]['depthlim'],
                                   marker=cfginfo[proc]['markertype'])
            elif cfginfo[proc]['action'] == 'markregion':
                basesmc.markRegion(cfginfo[proc]['extents'][0],
                                   cfginfo[proc]['extents'][1],
                                   cfginfo[proc]['extents'][2],
                                   cfginfo[proc]['extents'][3],
                                   marker=cfginfo[proc]['markertype'],
                                   depthlim=cfginfo[proc]['depthlim'])
            elif cfginfo[proc]['action'] == 'unmark':
                basesmc.unmarkCells(marker=cfginfo[proc]['markertype'],
                                   box=[cfginfo[proc]['extents'][0],
                                        cfginfo[proc]['extents'][1],
                                        cfginfo[proc]['extents'][2],
                                        cfginfo[proc]['extents'][3]],
                                   osbox=cfginfo[proc]['osbox'],
                                   thruzero=cfginfo[proc]['thruzero'])
                # deldry option removes dry cells as part of unmark action
                if cfginfo[proc]['deldry']:
                    basesmc.delCells(celltype='dry')
            if lp == len(cfginfo['marklist'])-1:
                # write the marked grid after final marking action
                basesmc.label = cfginfo[proc]['label']
                markfile = basesmc.writeNC(writedir=cfginfo[proc]['workdir'])
                # visualise the new marked grid
                grd.plotGridsmc(basesmc)

    elif cfginfo['action'] == 'tiergen':
        # load the parent grid
        print('')
        print('*** Loading SMC base grid')
        basefile = cfginfo['workdir'] + '/' + cfginfo['basefile']
        basesmcrm = grd.loadNCsmc(basefile)
        # generate the new tier grid
        print('')
        print('*** Creating tier for SMC grid')
        tiersmc = grd.createTiersmc(cfginfo['bathyfile'],
                                    basesmcrm,
                                    label=cfginfo['label'],
                                    mindepth=cfginfo['mindepth'],
                                    dlim=cfginfo['drydepthlim'], 
                                    drymin=cfginfo['drypcmin'], 
                                    drymax=cfginfo['drypcmax'], 
                                    bathytype=cfginfo['bathytype'],
                                    getpland=cfginfo['getpcland'],
                                    setadj=cfginfo['setadj'],
                                    deldry=cfginfo['deldry'])
        # save the new tier grid
        tierfile = tiersmc.writeNC(writedir=cfginfo['workdir'])
        # visualise the new tier grid
        grd.plotGridsmc(tiersmc)

    elif cfginfo['action'] == 'tiercombine':
        # load the parent grid
        print('')
        print('*** Loading SMC base grid')
        basefile = cfginfo['workdir'] + '/' + cfginfo['basefile']
        basesmcrm = grd.loadNCsmc(basefile)
        # load the tier
        print('')
        print('*** Loading SMC tier')
        tierfile = cfginfo['workdir'] + '/' + cfginfo['tierfile']
        tiersmc = grd.loadNCsmc(tierfile)
        # combine the grids
        print('')
        print('*** Combining tier and base')
        combbaset = grd.joinTiersmc(basesmcrm,
                                    tiersmc,
                                    cfginfo['bathyfile'],
                                    tiernext=cfginfo['tiernext'])
        combbaset.label = cfginfo['label']
        # save the combined grid
        combbasetfile = combbaset.writeNC(writedir=cfginfo['workdir'])
        # visualise the grid
        grd.plotGridsmc(combbaset, latlon=False)

    elif cfginfo['action'].lower() == 'writeww3':
        # load the grid for writing to ww3 format
        print('')
        print('*** Loading combined grid')
        gridfile = cfginfo['workdir'] + '/' + cfginfo['gridfile']
        combsmc = grd.loadNCsmc(gridfile)
        # remove any dry cells from the grid
        print('')
        print('*** Removing dry cells')
        combsmc.delCells(celltype='alldry')
        # visualize the grid
        grd.plotGridsmc(combsmc, latlon=False)
        # write out ww3 format files
        print('')
        print('*** Writing smc cells to text file')
        combsmc.sortCells()
        combsmc.writeWW3(writedir=cfginfo['writedir'],
                         mindepth=cfginfo['mindepth'],
                         writemindepth=cfginfo['writemindepth'],
                         arctic=cfginfo['arctic'],
                         arclat=cfginfo['arclat'])

    return


## main - read in the required action and run scripts accordingly
procname = sys.argv[1]
if len(sys.argv) > 2:
    cfgfile = sys.argv[2]
else:
    cfgfile = './gridgen.cfg'

cfginfo = readConfig(procname, cfgfile)
runGridgen(cfginfo)

