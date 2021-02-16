#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import modeller
import ppdg
import logging
log = logging.getLogger(__name__)

        !!!!!!!!!!!!!!!!!!
            Do not use
        !!!!!!!!!!!!!!!!!!

def ga341(wrkdir):
    """
        Get GA341 score for protein-protein. PDB must contain only two chains.
    """
    time_start = timer()
    log.info("Getting GA341 scoring...")
    cpx = os.path.join(wrkdir, 'complex.pdb')
    rec = os.path.join(wrkdir, 'receptor.pdb')
    lig = os.path.join(wrkdir, 'ligand.pdb')

    with open("ga341.out", "w") as fp:
        _stdout = sys.stdout
        sys.stdout = fp
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        cpx = modeller.scripts.complete_pdb(env, cpx)
        rec = modeller.scripts.complete_pdb(env, rec)
        lig = modeller.scripts.complete_pdb(env, lig)
        cpx.seq_id = 0.0
        rec.seq_id = 0.0
        lig.seq_id = 0.0
        cpx = cpx.assess_ga341()
        rec = rec.assess_ga341()
        lig = lig.assess_ga341()
        print([ c-r-l for c,r,l in zip(cpx, rec, lig)])
        quit()
        score = cpx.assess_ga341() - rec.assess_ga341() - lig.assess_ga341()
        sys.stdout.flush()
        sys.stdout = _stdout

    desc = dict()
    desc['GA341'] = score
    time_end = timer()
    desc['>TIME_GA341'] = time_end - time_start
    return desc

if __name__=='__main__':
    print(ga341(sys.argv[1]))
