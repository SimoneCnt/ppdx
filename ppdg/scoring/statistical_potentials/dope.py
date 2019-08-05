#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import modeller
import ppdg
import logging
log = logging.getLogger(__name__)

def dope(wrkdir):
    """
        Get DOPE and DOPE-HRS scores for protein-protein.
        PDB must contain only two chains.
        [1] M. Shen and A. Sali, "Statistical potential for assessment and 
            prediction of protein structures", Protein Science, vol. 15, 
            no. 11, pp. 2507-2524, 2006.
    """
    time_start = timer()
    log.info("Getting DOPE and DOPE-HR scoring...")
    cpx = os.path.join(wrkdir, 'complexAB.pdb')

    with open("dope.out", "w") as fp:
        _stdout = sys.stdout
        sys.stdout = fp
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, cpx)
        cpx = modeller.selection(mdl.chains[:])
        lig = modeller.selection(mdl.chains[0])
        rec = modeller.selection(mdl.chains[1])
        dope = cpx.assess_dope() - rec.assess_dope() - lig.assess_dope()
        dopehr = cpx.assess_dopehr() - rec.assess_dopehr() - lig.assess_dopehr()
        sys.stdout.flush()
        sys.stdout = _stdout

    desc = dict()
    desc['DOPE'] = dope
    desc['DOPE-HR'] = dopehr
    time_end = timer()
    desc['>TIME_DOPE'] = time_end - time_start
    return desc

