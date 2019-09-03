#!/usr/bin/env python3

import os, sys
import numpy as np
from timeit import default_timer as timer
import ppdg
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
    pdb = ppdg.Pdb('ligand.pdb')
    pdb.remove_hydrogens()
    pdb.write('ligand_firedock_noh.pdb')
    pdb = ppdg.Pdb('receptor.pdb')
    pdb.remove_hydrogens()
    pdb.write('receptor_firedock_noh.pdb')

    # Reduce
    fd_reduce = "%s/PDBPreliminaries/reduce.2.21.030604 -DB %s/PDBPreliminaries/reduce_het_dict.txt -OH -HIS -NOADJust -NOROTMET" % (ppdg.FIREDOCK, ppdg.FIREDOCK)
    stdout, stderr, ret = ppdg.tools.execute(fd_reduce+' receptor_firedock_noh.pdb')
    with open('firedock_reduce_receptor.out', 'w') as fp:
        fp.write(stderr)
    with open('firedock_receptor.pdb', 'w') as fp:
        fp.write(stdout)
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock reduce receptor failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))
    stdout, stderr, ret = ppdg.tools.execute(fd_reduce+' ligand_firedock_noh.pdb')
    with open('firedock_reduce_ligand.out', 'w') as fp:
        fp.write(stderr)
    with open('firedock_ligand.pdb', 'w') as fp:
        fp.write(stdout)
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock reduce ligand failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))

    # Score
    with open('firedock.trans', 'w') as fp:
        fp.write("1 0.0 0.0 0.0 0.0 0.0 0.0")
    stdout, stderr, ret = ppdg.tools.execute("%s/buildFireDockParams.pl firedock_receptor.pdb firedock_ligand.pdb U U Default firedock.trans firedock_build.out 1 50 0.85 0 firedock_parameters.dat" % ppdg.FIREDOCK)
    with open("firedock_build.err", 'w') as fp:
        fp.write(stdout+'\n'+stderr)
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock build param failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))
    stdout, stderr, ret = ppdg.tools.execute("%s/runFireDock.pl firedock_parameters.dat" % (ppdg.FIREDOCK))
    with open('firedock.log', 'w') as fp:
        fp.write(stdout+'\n'+stderr)
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("FireDock scoring failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))

    # Get values
    with open('firedock_build.out.ref', 'r') as fp:
        data = fp.readlines()[-1].split('|')
    if data[5].strip()=='glob':
        raise ValueError('Problem with FireDock in %s' % (wrkdir))
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
    ppdg.readconfig('config-ppdg.ini')
    print(firedock(sys.argv[1]))

