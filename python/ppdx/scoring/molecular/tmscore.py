#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import ppdx
import logging
log = logging.getLogger(__name__)

def tmscore(wrkdir):
    """
        TM-align score - structural similarity between the model and the template.
        Useful as a diagnostic test tomake sure the models are generated correctly.
    """
    time_start = timer()
    log.info("Getting TM-score...")
    if not os.path.isfile(os.path.join(wrkdir, 'complex.pdb')):
        raise ValueError('File complex.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'template.pdb')):
        raise ValueError('File template.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)
    ret = ppdx.tools.execute("%s %s %s -outfmt 2 >tmscore.out 2>&1" % (ppdx.TMALIGN, 'template.pdb', 'complex.pdb'))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("TM-align failed!")
    with open('tmscore.out') as fp:
        splt = fp.readlines()[-2].split()
        score = min(float(splt[2]), float(splt[3]))
    os.chdir(basepath)
    time_end = timer()
    desc = dict()
    desc['TMscore'] = score
    desc['>TIME_TMscore'] = time_end - time_start
    return desc

