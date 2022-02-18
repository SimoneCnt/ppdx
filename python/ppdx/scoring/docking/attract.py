#!/usr/bin/env python3

import os, sys
import numpy as np
from timeit import default_timer as timer
import ppdx
import logging
log = logging.getLogger(__name__)

def attract(wrkdir):
    """
        Calculate ATTRACT scoring.
    """
    time_start = timer()
    log.info("Getting ATTRACT scoring...")
    basepath = os.getcwd()
    os.chdir(wrkdir)

    # Create input file
    with open('attract.dat', 'w') as fp:
        fp.write("## SimoRun\n")
        fp.write("#pivot 1 0. 0. 0. \n")
        fp.write("#pivot 2 0. 0. 0. \n")
        fp.write("#centered receptor: false \n")
        fp.write("#centered ligands: false \n")
        fp.write("#1\n")
        fp.write("   0.000000 0.000000 0.000000 0.0000 0.0000 0.0000\n")
        fp.write("   0.000000 0.000000 0.000000 0.0000 0.0000 0.0000\n")

    # List of commands to execute
    cmds = [
        "%s/bin/shm-clean >/dev/null 2>&1" % (ppdx.ATTRACT),
        "python2 %s/allatom/aareduce.py receptor.pdb receptor-aa.pdb --chain A --pdb2pqr >attract_receptor_reduce.out 2>&1" % (ppdx.ATTRACT),
        "python2 %s/tools/reduce.py receptor-aa.pdb receptorr.pdb --chain A >attract_receptor_reduce2.out 2>&1" % (ppdx.ATTRACT),
        "python2 %s/allatom/aareduce.py ligand.pdb ligand-aa.pdb --chain B --pdb2pqr >attract_ligand_reduce.out 2>&1" % (ppdx.ATTRACT),
        "python2 %s/tools/reduce.py ligand-aa.pdb ligandr.pdb --chain B >attract_ligand_reduce2.out" % (ppdx.ATTRACT),
        "%s/bin/attract attract.dat %s/attract.par receptorr.pdb ligandr.pdb --score --fix-receptor --rcut 50.0 >attract.score 2>&1" % (ppdx.ATTRACT, ppdx.ATTRACT),
        "%s/bin/shm-clean >/dev/null 2>&1" % (ppdx.ATTRACT)
    ]

    # Exec everything...
    for cmd in cmds:
        ret = ppdx.tools.execute(cmd)
        if ret!=0:
            os.chdir(basepath)
            raise ValueError("ATTRACT command failed!\n<%s>" % (cmd))

    # Get the score
    with open("attract.score", 'r') as fp:
        score = float(fp.readline().split()[1])
    os.chdir(basepath)
    time_end = timer()
    desc = {'ATTRACT':score, '>TIME_ATTRACT':time_end-time_start}
    return desc

if __name__=='__main__':
    ppdx.config.cread('config-ppdx.ini')
    print(attract(sys.argv[1]))

