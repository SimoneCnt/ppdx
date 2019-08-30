#!/usr/bin/env python3

import os, sys
import numpy as np
from timeit import default_timer as timer
import ppdg
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
        ("%s/bin/shm-clean" % (ppdg.ATTRACT), None),
        ("python2 %s/allatom/aareduce.py receptor.pdb receptor-aa.pdb --chain A --pdb2pqr" % (ppdg.ATTRACT), "attract_receptor_reduce.out"),
        ("python2 %s/tools/reduce.py receptor-aa.pdb receptorr.pdb --chain A" % (ppdg.ATTRACT), "attract_receptor_reduce2.out"),
        ("python2 %s/allatom/aareduce.py ligand.pdb ligand-aa.pdb --chain B --pdb2pqr" % (ppdg.ATTRACT), "attract_ligand_reduce.out"),
        ("python2 %s/tools/reduce.py ligand-aa.pdb ligandr.pdb --chain B" % (ppdg.ATTRACT), "attract_ligand_reduce2.out"),
        ("%s/bin/attract attract.dat %s/attract.par receptorr.pdb ligandr.pdb --score --fix-receptor --rcut 50.0" % (ppdg.ATTRACT, ppdg.ATTRACT), "attract.score"),
        ("%s/bin/shm-clean" % (ppdg.ATTRACT), None)
    ]

    # Exec everything...
    for cmd, outfile in cmds:
        stdout, stderr, ret = ppdg.tools.execute(cmd)
        if ret!=0:
            os.chdir(basepath)
            raise ValueError("ATTRACT command failed!\n<%s>\nReturned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (cmd, ret, stdout, stderr))
        if outfile:
            with open(outfile, 'w') as fp:
                fp.write(stdout+'\n'+stderr)

    # Get the score
    with open("attract.score", 'r') as fp:
        score = float(fp.readline().split()[1])
    os.chdir(basepath)
    time_end = timer()
    desc = {'ATTRACT':score, '>TIME_ATTRACT':time_end-time_start}
    return desc

if __name__=='__main__':
    ppdg.readconfig('config-ppdg.ini')
    print(attract(sys.argv[1]))

