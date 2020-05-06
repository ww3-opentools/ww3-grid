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
# run_gebco_reduce.py
#
# PURPOSE:
#  Runs gebco_reduce.py functions based on user defined configuration file inputs
#
# REVISION HISTORY:
#
# A. Saulter; Met Office; May-2020 Version:1.0
#  Code prepared for initial release on github
#
#==================================================================================

import matplotlib.pyplot as plt
import numpy as np
import configparser as cfg
import sys

from gebco_reduce import reduceGEBCO
from gebco_reduce import interpGEBCO
from gebco_reduce import correctLakesBathyfromfile
from gebco_reduce import plotGridfromfile

def readConfig(action, cfgfile):
    '''Reads in relevant content of config file according to action'''

    if action == 'reduce':
        # set defaults
        cfginfo = {'scalefac':6,
                   'depthmin':0.0,
                   'region':None,
                   'extents':None,
                   'pltchk':True,
                   'correctlakes':False,
                   'gebcofile':'GEBCO_2014_2D.nc',
                   'datadir':'.'}
        # read updated info from the configuration file
        print('[INFO] Searching for configuration information from: %s' %cfgfile)
        config = cfg.RawConfigParser()
        try:
            config.read(cfgfile)
            if config.has_option(action,'scalefac'):
                cfginfo['scalefac'] = np.int(config.get(action,'scalefac'))
            if config.has_option(action,'depthmin'):
                cfginfo['depthmin'] = np.float(config.get(action,'depthmin'))
            if config.has_option(action,'region'):
                if config.get(action,'region').lower() != 'none':
                    cfginfo['region'] = config.get(action,'region')
            if config.has_option(action,'extents'):
                if config.get(action,'extents').lower() != 'none':
                    extentstr = config.get(action,'extents').split(',')
                    cfginfo['extents'] = [np.float(extentstr[0]), np.float(extentstr[1]),
                                          np.float(extentstr[2]), np.float(extentstr[3])]
            if config.has_option(action,'pltchk'):
                if config.get(action,'pltchk').lower() == 'false':
                    cfginfo['pltchk'] = False
            if config.has_option(action,'correctlakes'):
                if config.get(action,'correctlakes').lower() == 'true':
                    cfginfo['correctlakes'] = True
            if config.has_option(action,'gebcofile'):
                cfginfo['gebcofile'] = config.get(action,'gebcofile')
            if config.has_option(action,'datadir'):
                cfginfo['datadir'] = config.get(action,'datadir')
        except:
            print('[WARN] Configuration file %s missing; returning defaults' %cfgfile)

    if action == 'interpolate':
        # set defaults
        cfginfo = {'dx':0.5,
                   'dy':0.5,
                   'depthmin':0.0,
                   'extents':None,
                   'pltchk':True,
                   'correctlakes':False,
                   'gebcofile':'GEBCO_2014_2D.nc',
                   'datadir':'.'}
        # read updated info from the configuration file
        print('[INFO] Searching for configuration information from: %s' %cfgfile)
        config = cfg.RawConfigParser()
        try:
            config.read(cfgfile)
            if config.has_option(action,'dx'):
                cfginfo['dx'] = np.float(config.get(action,'dx'))
            if config.has_option(action,'dy'):
                cfginfo['dy'] = np.float(config.get(action,'dy'))
            if config.has_option(action,'depthmin'):
                cfginfo['depthmin'] = np.float(config.get(action,'depthmin'))
            if config.has_option(action,'region'):
                if config.get(action,'region').lower() != 'none':
                    cfginfo['region'] = config.get(action,'region')
            if config.has_option(action,'extents'):
                if config.get(action,'extents').lower() != 'none':
                    extentstr = config.get(action,'extents').split(',')
                    cfginfo['extents'] = [np.float(extentstr[0]), np.float(extentstr[1]),
                                          np.float(extentstr[2]), np.float(extentstr[3])]
            if config.has_option(action,'pltchk'):
                if config.get(action,'pltchk').lower() == 'false':
                    cfginfo['pltchk'] = False
            if config.has_option(action,'correctlakes'):
                if config.get(action,'correctlakes').lower() == 'true':
                    cfginfo['correctlakes'] = True
            if config.has_option(action,'gebcofile'):
                cfginfo['gebcofile'] = config.get(action,'gebcofile')
            if config.has_option(action,'datadir'):
                cfginfo['datadir'] = config.get(action,'datadir')
        except:
            print('[WARN] Configuration file %s missing; returning defaults' %cfgfile)

    elif action == 'correct':
        # set defaults
        cfginfo = {'depthmin':0.0,
                   'pltchk':True,
                   'ncfile':'GEBCO_2014_2D_reduced.nc',
                   'datadir':'.',
                   'removesmall':None,
                   'caspianonly':False}
        # read updated info from the configuration file
        print('[INFO] Searching for configuration information from: %s' %cfgfile)
        config = cfg.RawConfigParser()
        try:
            config.read(cfgfile)
            if config.has_option(action,'depthmin'):
                cfginfo['depthmin'] = np.float(config.get(action,'depthmin'))
            if config.has_option(action,'pltchk'):
                if config.get(action,'pltchk').lower() == 'false':
                    cfginfo['pltchk'] = False
            if config.has_option(action,'ncfile'):
                cfginfo['ncfile'] = config.get(action,'ncfile')
            if config.has_option(action,'datadir'):
                cfginfo['datadir'] = config.get(action,'datadir')
            if config.has_option(action,'removesmall'):
                if config.get(action,'removesmall').lower() != 'none':
                    cfginfo['removesmall'] = np.int(config.get(action,'removesmall'))
            if config.has_option(action,'caspianonly'):
                if config.get(action,'caspianonly').lower() == 'true':
                    cfginfo['caspianonly'] = True
        except:
            print('[WARN] Configuration file %s missing; returning defaults' %cfgfile)

    elif action == 'plot':
        # set defaults
        cfginfo = {'depthmin':5.0,
                   'depthmax':-500.0,
                   'ncfile':'GEBCO_2014_2D_reduced.nc',
                   'datadir':'.',
                   'usedepths':True,
                   'uselandsea':True}
        # read updated info from the configuration file
        print('[INFO] Searching for configuration information from: %s' %cfgfile)
        config = cfg.RawConfigParser()
        try:
            config.read(cfgfile)
            if config.has_option(action,'depthmin'):
                cfginfo['depthmin'] = np.float(config.get(action,'depthmin'))
            if config.has_option(action,'depthmax'):
                cfginfo['depthmax'] = np.float(config.get(action,'depthmax'))
            if config.has_option(action,'ncfile'):
                cfginfo['ncfile'] = config.get(action,'ncfile')
            if config.has_option(action,'datadir'):
                cfginfo['datadir'] = config.get(action,'datadir')
            if config.has_option(action,'usedepths'):
                if config.get(action,'usedepths').lower() == 'false':
                    cfginfo['usedepths'] = False
            if config.has_option(action,'uselandsea'):
                if config.get(action,'uselandsea').lower() == 'false':
                    cfginfo['uselandsea'] = False
        except:
            print('[WARN] Configuration file %s missing; returning defaults' %cfgfile)

    return cfginfo        


## main - read in the required action and run scripts accordingly
action  = sys.argv[1]
if len(sys.argv) > 2:
    cfgfile = sys.argv[2]
else:
    cfgfile = './gebco_reduce.cfg'

if action.lower() == 'reduce':
    cfginfo = readConfig('reduce', cfgfile)
    reduceGEBCO(scalefac=cfginfo['scalefac'], 
                 depthmin=cfginfo['depthmin'],
                 region=cfginfo['region'],
                 cutout=cfginfo['extents'],
                 pltchk=cfginfo['pltchk'],
                 correctlakes=cfginfo['correctlakes'],
                 gebcofile=cfginfo['gebcofile'],
                 datadir=cfginfo['datadir'])

elif action.lower() == 'interpolate':
    cfginfo = readConfig('interpolate', cfgfile)
    interpGEBCO(dx=cfginfo['dx'], 
                 dy=cfginfo['dy'], 
                 depthmin=cfginfo['depthmin'],
                 region=cfginfo['region'],
                 cutout=cfginfo['extents'],
                 pltchk=cfginfo['pltchk'],
                 correctlakes=cfginfo['correctlakes'],
                 gebcofile=cfginfo['gebcofile'],
                 datadir=cfginfo['datadir'])

elif action.lower() == 'correct':
    cfginfo = readConfig('correct', cfgfile)
    correctLakesBathyfromfile(cfginfo['ncfile'], 
                               datadir=cfginfo['datadir'],
                               depthmin=cfginfo['depthmin'],
                               removesmall=cfginfo['removesmall'],
                               caspianonly=cfginfo['caspianonly'],
                               pltchk=cfginfo['pltchk'])

elif action.lower() == 'plot':
    cfginfo = readConfig('plot', cfgfile)
    plotGridfromfile(cfginfo['ncfile'], 
                  datadir=cfginfo['datadir'],
                  usedepths=cfginfo['usedepths'],
                  uselandsea=cfginfo['uselandsea'],
                  depthmin=cfginfo['depthmin'],
                  depthmax=cfginfo['depthmax'])
