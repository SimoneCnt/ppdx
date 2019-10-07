#!/usr/bin/env python3

import os
import ppdg
import logging
log = logging.getLogger(__name__)

def lisa(wrkdir, numthreads=None):
    """
        Calculate LISA scoring.
    """
    if not os.path.isfile(os.path.join(wrkdir, 'ligand.pdb')):
        raise ValueError('File ligand.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'receptor.pdb')):
        raise ValueError('File receptor.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)

    # Getting the number fo cores to use
    if not numthreads:
        numthreads = os.cpu_count()
        if not numthreads:
            numthreads = 1

    log.info("Getting LISA scoring using %d threads..." % (numthreads))

    ret = ppdg.tools.execute('LISA-1.0.py receptor.pdb ligand.pdb %d >lisa.out 2>&1' % (numthreads))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("LISA failed!")

    outfile = os.path.join(wrkdir, 'LISA-output', 'ligand', 'LISA_bindingAff.txt')
    if not os.path.isfile(outfile):
        raise ValueError('LISA output file %s not found!' % (outfile))
    with open(outfile) as fp:
        fp.readline()
        score = float(fp.readline().split()[1])

    os.chdir(basepath)
    return {'LISA' : score}

