#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def zrank(wrkdir):
    """
        Calculate ZRANK scoring.
    """
    time_start = timer()
    log.info("Getting ZRANK scoring...")
    if not os.path.isfile(os.path.join(wrkdir, 'complexAB.pdb')):
        raise ValueError('File complexAB.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)
    with open("zrank", 'w') as fp:
        fp.write('complexAB.pdb')
    stdout, stderr, ret = ppdg.tools.execute("%s zrank" % (os.path.join(ppdg.ZRANK, 'zrank')))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("ZRANK failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))
    with open("zrank.zr.out", 'r') as fp:
        zrank = float(fp.readline().split()[1])
    os.chdir(basepath)
    desc = {'ZRANK' : zrank}
    time_end = timer()
    desc['>TIME_ZRANK'] = time_end - time_start
    return desc

def zrank2(wrkdir):
    """
        Calculate ZRANK2 scoring.
    """
    time_start = timer()
    log.info("Getting ZRANK2 scoring...")
    if not os.path.isfile(os.path.join(wrkdir, 'complexAB.pdb')):
        raise ValueError('File complexAB.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)
    with open("zrank2", 'w') as fp:
        fp.write('complexAB.pdb')
    stdout, stderr, ret = ppdg.tools.execute("%s -R zrank2" % (os.path.join(ppdg.ZRANK, 'zrank')))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("ZRANK failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))
    with open("zrank2.zr.out", 'r') as fp:
        zrank2 = float(fp.readline().split()[1])
    os.chdir(basepath)
    desc = {'ZRANK2' : zrank2}
    time_end = timer()
    desc['>TIME_ZRANK2'] = time_end - time_start
    return desc

