#!/usr/bin/env python3

import os
import configparser
import ppdg

def readconfig(fname):
    '''
        Set some needed variables.
    '''
    defaults = dict(
            WRKDIR      = os.path.join(os.getcwd(), 'models'),
            ZRANK       = "",
            RFSPP       = "",
            FOLDX       = "",
            IPOT        = "",
            PYDOCK      = "",
            CHARMM      = "",
            ATTRACT     = "",
            PDBDIR      = "",
            FFPATH      = "",
            ROSETTA     = "",
            ROSETTABIN  = ""
        )
    config = configparser.ConfigParser(defaults)
    config.read(fname)
    ppdg.WRKDIR = config.get('ppdg', 'WRKDIR')
    ppdg.ZRANK  = config.get('ppdg', 'ZRANK')
    ppdg.RFSPP  = config.get('ppdg', 'RFSPP')
    ppdg.FOLDX  = config.get('ppdg', 'FOLDX')
    ppdg.IPOT   = config.get('ppdg', 'IPOT')
    ppdg.PYDOCK = config.get('ppdg', 'PYDOCK')
    ppdg.CHARMM = config.get('ppdg', 'CHARMM')
    ppdg.ATTRACT= config.get('ppdg', 'ATTRACT')
    ppdg.FIREDOCK= config.get('ppdg', 'FIREDOCK')
    ppdg.PDBDIR = config.get('ppdg', 'PDBDIR')
    ppdg.FFPATH = config.get('ppdg', 'FFPATH')
    ppdg.ROSETTA = config.get('ppdg', 'ROSETTA')
    ppdg.ROSETTABIN = config.get('ppdg', 'ROSETTABIN')

def printconfig():
    '''
        Print current settings.
    '''
    for key, value in ppdg.__dict__.items():
        if not (key.startswith('__') or key.startswith('_')) and key==key.upper():
            print("%-10s = %s" % (key, value))

