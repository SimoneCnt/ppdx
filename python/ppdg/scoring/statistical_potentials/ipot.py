#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def ipot_core(wrkdir, pote):
    """
        Calculate iPot scoring [1].
        Code from https://github.com/gjoni/iPot
        Commit 03ffb173dbca25a31bc3a8f3ebe0ce9c960578e5 )
        [1] I. Anishchenko, P. J. Kundrotas, and I. A. Vakser, "Contact 
            Potential for Structure Prediction of Proteins and Protein 
            Complexes from Potts Model", Biophysical Journal, vol. 115, 
            no. 5, pp. 809-821, 2018.
    """
    time_start = timer()
    log.info("Getting iPot %s scoring..." % (pote))
    if not os.path.isfile(os.path.join(wrkdir, 'complex.pdb')):
        raise ValueError('File complex.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'ligand.pdb')):
        raise ValueError('File ligand.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'receptor.pdb')):
        raise ValueError('File receptor.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)
    outfile = 'ipot_%s.out' % (pote)
    ret = ppdg.tools.execute("%s -r receptor.pdb -l ligand.pdb >%s 2>&1" % (os.path.join(ppdg.IPOT, pote), outfile))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("iPot with potential %s failed!" % (pote))
    with open(outfile) as fp:
        dG = float(fp.readlines()[-1].split()[1])
    os.chdir(basepath)
    time_end = timer()
    desc = dict()
    desc['ipot_'+pote] = dG
    desc['>TIME_ipot_'+pote] = time_end - time_start
    return desc

def ipot_aace167(wrkdir):
    return ipot_core(wrkdir, 'aace167')

def ipot_aace18(wrkdir):
    return ipot_core(wrkdir, 'aace18')

def ipot_aace20(wrkdir):
    return ipot_core(wrkdir, 'aace20')

def ipot_rrce20(wrkdir):
    return ipot_core(wrkdir, 'rrce20')


