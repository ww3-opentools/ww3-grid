# run gridgen actions based on user inputs

import matplotlib.pyplot as plt
import numpy as np
import configparser as cfg
import sys

import gridgen as grd

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
                   'writemindepth':False}
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

    return cfginfo        


## main - read in the required action and run scripts accordingly
procname = sys.argv[1]
if len(sys.argv) > 2:
    cfgfile = sys.argv[2]
else:
    cfgfile = './gridgen.cfg'

cfginfo = readConfig(procname, cfgfile)

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
    # write staging file
    basefile = basesmc.writeNC(writedir=cfginfo['workdir'])
    print('Data written to %s' %basefile)
    # visualize the grid
    grd.plotGridsmc(basesmc, latlon=False)

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
    print('Data written to %s' %tierfile)
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
    # save the grid
    combbasetfile = combbaset.writeNC(writedir=cfginfo['workdir'])
    print('Data written to %s' %combbasetfile)
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
    combsmc.delDryCells(celltype='alldry')
    # visualize the grid
    grd.plotGridsmc(combsmc, latlon=False)
    # write out ww3 format files
    print('')
    print('*** Writing smc cells to text file')
    combsmc.sortCells()
    combsmc.writeWW3(writedir=cfginfo['writedir'],
                     mindepth=cfginfo['mindepth'],
                     writemindepth=cfginfo['writemindepth'])
