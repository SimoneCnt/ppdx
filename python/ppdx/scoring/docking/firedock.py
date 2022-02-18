#!/usr/bin/env python3

import os, sys
import numpy as np
from timeit import default_timer as timer
import ppdx
import logging
log = logging.getLogger(__name__)

def firedock(wrkdir):
    """
        Calculate FireDock scoring.
    """
    time_start = timer()
    log.info("Getting FireDock scoring...")
    basepath = os.getcwd()
    os.chdir(wrkdir)

    # Remove hydrogens
    pdb = ppdx.Pdb('ligand.pdb')
    pdb.remove_hydrogens()
    pdb.write('ligand_firedock_noh.pdb')
    pdb = ppdx.Pdb('receptor.pdb')
    pdb.remove_hydrogens()
    pdb.write('receptor_firedock_noh.pdb')

    # Reduce
    fd_reduce = "%s/PDBPreliminaries/reduce.2.21.030604 -DB %s/PDBPreliminaries/reduce_het_dict.txt -OH -HIS -NOADJust -NOROTMET" % (ppdx.FIREDOCK, ppdx.FIREDOCK)
    ret = ppdx.tools.execute(fd_reduce+' receptor_firedock_noh.pdb >firedock_receptor.pdb 2>firedock_reduce_receptor.out')
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock reduce receptor failed!")
    ret = ppdx.tools.execute(fd_reduce+' ligand_firedock_noh.pdb >firedock_ligand.pdb 2>firedock_reduce_ligand.out')
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock reduce ligand failed!")

    # Score
    with open('firedock.trans', 'w') as fp:
        fp.write("1 0.0 0.0 0.0 0.0 0.0 0.0")
    ret = ppdx.tools.execute("%s/buildFireDockParams.pl firedock_receptor.pdb firedock_ligand.pdb U U Default firedock.trans firedock_build.out 1 50 0.85 0 firedock_parameters.dat >firedock_build.err 2>&1" % ppdx.FIREDOCK)
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock build param failed!")
    ret = ppdx.tools.execute("%s/runFireDock.pl firedock_parameters.dat >firedock.log 2>&1" % (ppdx.FIREDOCK))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock scoring failed!")

    # Get values
    with open('firedock_build.out.ref', 'r') as fp:
        data = fp.readlines()[-1].split('|')
    if data[5].strip()=='glob':
        log.warning('Problem with FireDock in %s' % (wrkdir))
        data = [float('nan')]*20
    desc = dict()
    desc['FireDock']          = float(data[5])
    desc['FireDock_aVdW']     = float(data[6])
    desc['FireDock_rVdW']     = float(data[7])
    desc['FireDock_ACE']      = float(data[8])
    desc['FireDock_inside']   = float(data[9])
    desc['FireDock_aElec']    = float(data[10])
    desc['FireDock_rElec']    = float(data[11])
    desc['FireDock_laElec']   = float(data[12])
    desc['FireDock_lrElec']   = float(data[13])
    desc['FireDock_hb']       = float(data[14])
    desc['FireDock_piS']      = float(data[15])
    desc['FireDock_catpiS']   = float(data[16])
    desc['FireDock_aliph']    = float(data[17])
    os.chdir(basepath)
    time_end = timer()
    desc['>TIME_FireDock'] = time_end - time_start
    return desc

if __name__=='__main__':
    ppdx.config.cread('config-ppdx.ini')
    print(firedock(sys.argv[1]))

