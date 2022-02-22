#!/usr/bin/env python3

import os
import configparser
import ppdx

def cread(fname):
    defaults = dict(
            WRKDIR      = os.path.join(os.getcwd(), 'models'),
            ZRANK       = "",
            RFSPP       = "",
            FOLDX       = "",
            IPOT        = "",
            PYDOCK      = "",
            CHARMM      = "",
            ATTRACT     = "",
            FIREDOCK    = "",
            PDBDIR      = "",
            FFPATH      = "",
            ROSETTA     = "",
            ROSETTABIN  = ""
        )
    config = configparser.ConfigParser(defaults)
    config.read(fname)
    ppdx.WRKDIR = config.get('ppdx', 'WRKDIR')
    ppdx.ZRANK  = config.get('ppdx', 'ZRANK')
    ppdx.RFSPP  = config.get('ppdx', 'RFSPP')
    ppdx.FOLDX  = config.get('ppdx', 'FOLDX')
    ppdx.IPOT   = config.get('ppdx', 'IPOT')
    ppdx.PYDOCK = config.get('ppdx', 'PYDOCK')
    ppdx.CHARMM = config.get('ppdx', 'CHARMM')
    ppdx.ATTRACT= config.get('ppdx', 'ATTRACT')
    ppdx.FIREDOCK= config.get('ppdx', 'FIREDOCK')
    ppdx.PDBDIR = config.get('ppdx', 'PDBDIR')
    ppdx.FFPATH = config.get('ppdx', 'FFPATH')
    ppdx.ROSETTA = config.get('ppdx', 'ROSETTA')
    ppdx.ROSETTABIN = config.get('ppdx', 'ROSETTABIN')

def cget():
    config = dict()
    config['WRKDIR']     = ppdx.WRKDIR
    config['ZRANK']      = ppdx.ZRANK
    config['RFSPP']      = ppdx.RFSPP
    config['FOLDX']      = ppdx.FOLDX
    config['IPOT']       = ppdx.IPOT
    config['PYDOCK']     = ppdx.PYDOCK
    config['CHARMM']     = ppdx.CHARMM
    config['ATTRACT']    = ppdx.ATTRACT
    config['FIREDOCK']   = ppdx.FIREDOCK
    config['PDBDIR']     = ppdx.PDBDIR
    config['FFPATH']     = ppdx.FFPATH
    config['ROSETTA']    = ppdx.ROSETTA
    config['ROSETTABIN'] = ppdx.ROSETTABIN
    return config

def cset(config):
    ppdx.WRKDIR     = config['WRKDIR']
    ppdx.ZRANK      = config['ZRANK']
    ppdx.RFSPP      = config['RFSPP']
    ppdx.FOLDX      = config['FOLDX']
    ppdx.IPOT       = config['IPOT']
    ppdx.PYDOCK     = config['PYDOCK']
    ppdx.CHARMM     = config['CHARMM']
    ppdx.ATTRACT    = config['ATTRACT']
    ppdx.FIREDOCK   = config['FIREDOCK']
    ppdx.PDBDIR     = config['PDBDIR']
    ppdx.FFPATH     = config['FFPATH']
    ppdx.ROSETTA    = config['ROSETTA']
    ppdx.ROSETTABIN = config['ROSETTABIN']
    return

def cprint():
    '''
        Print current settings.
    '''
    c = cget()
    for key, value in c.items():
        print("%-10s = %s" % (key, value))

