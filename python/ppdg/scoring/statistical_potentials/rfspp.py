#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def rfspp_core_core(pdb, pote):
    """
        Compute the binding affinity using RF_HA_SRS and RF_CB_SRS_OD scoring functions.
        [1] D. Rykunov and A. Fiser, "Effects of amino acid composition, finite 
            size of proteins, and sparse statistics on distance-dependent 
            statistical pair potentials", Proteins: Structure, Function, and 
            Bioinformatics, vol. 67, no. 3, pp. 559-568, 2007.
        [2] D. Rykunov and A. Fiser, "New statistical potential for quality 
            assessment of protein models and a survey of energy functions", 
            BMC Bioinformatics, vol. 11, p. 128, 2010.
    """
    rfspp = os.path.join(ppdg.RFSPP, "calc_energy")
    outfile = 'rfspp_%s.out' % (pote)
    ret = ppdg.tools.execute("%s %s %s >%s 2>&1" % (rfspp, pdb, os.path.join(ppdg.RFSPP, pote), outfile))
    if ret!=0:
        raise ValueError("calc_energy with potential %s failed!" % (pote))
    with open(outfile) as fp:
        dg = float(fp.readlines()[-1])
    return dg

def rfspp_core(wrkdir, pote):
    """
        Calculate RFSPP scoring.
    """
    time_start = timer()
    log.info("Getting %s scoring..." % (pote))
    if not os.path.isfile(os.path.join(wrkdir, 'complex.pdb')):
        raise ValueError('File complex.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'ligand.pdb')):
        raise ValueError('File ligand.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'receptor.pdb')):
        raise ValueError('File receptor.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)
    cpx = rfspp_core_core('complex.pdb', pote)
    lig = rfspp_core_core('ligand.pdb', pote)
    rec = rfspp_core_core('receptor.pdb', pote)
    os.chdir(basepath)
    time_end = timer()
    return cpx-lig-rec, time_end-time_start

def rf_cb_srs_od(wrkdir):
    """
        Calculate RF_CB_SRS_OD scoring.
    """
    dg, time = rfspp_core(wrkdir, 'RF_CB_SRS_OD')
    desc = dict()
    desc['RF_CB_SRS_OD'] = dg
    desc['>TIME_RF_CB_SRS_OD'] = time
    return desc

def rf_ha_srs(wrkdir):
    """
        Calculate RF_HA_SRS scoring.
    """
    dg, time = rfspp_core(wrkdir, 'RF_HA_SRS')
    desc = dict()
    desc['RF_HA_SRS'] = dg
    desc['>TIME_RF_HA_SRS'] = time
    return desc

