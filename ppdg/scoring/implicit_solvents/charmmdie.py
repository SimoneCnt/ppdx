#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def run_die(fname, mode):
    if mode not in ['cdie', 'rdie']:
        raise ValueError('Mode must be either cdie or rdie')
    cmd = "%s sysname=%s -i %s.inp" % (os.path.join(ppdg.CHARMM, 'charmm'), fname, mode)
    out, err, ret = ppdg.tools.execute(cmd)
    with open('%s-%s.out' % (fname, mode), 'w') as fp:
        fp.write(out)
    with open('%s-%s.err' % (fname, mode), 'w') as fp:
        fp.write(err)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))
    for line in out.splitlines():
        if line.startswith('ENER EXTERN>'):
            _, _, vdw, elec, _, _, _ = line.split()
    return float(elec), float(vdw)

def cdie(wrkdir):
    """
        Get CDIE solvation term for protein-protein as implemented in CHARMM
    """
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting CDIE scoring...")
    ppdg.link_data('cdie.inp')

    cpx_elec, cpx_vdw = run_die('complex-chm', 'cdie')
    rec_elec, rec_vdw = run_die('receptor-chm', 'cdie')
    lig_elec, lig_vdw = run_die('ligand-chm', 'cdie')

    desc = dict()
    desc['CDIE_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['CDIE_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    time_end = timer()
    desc['>TIME_CDIE'] = time_end - time_start
    os.chdir(basepath)

    return desc


def rdie(wrkdir):
    """
        Get RDIE solvation term for protein-protein as implemented in CHARMM
    """
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting RDIE scoring...")
    ppdg.link_data('rdie.inp')

    cpx_elec, cpx_vdw = run_die('complex-chm', 'rdie')
    rec_elec, rec_vdw = run_die('receptor-chm', 'rdie')
    lig_elec, lig_vdw = run_die('ligand-chm', 'rdie')

    desc = dict()
    desc['RDIE_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['RDIE_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    time_end = timer()
    desc['>TIME_RDIE'] = time_end - time_start
    os.chdir(basepath)

    return desc


