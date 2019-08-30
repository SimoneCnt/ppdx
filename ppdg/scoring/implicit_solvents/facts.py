#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def run_facts(fname):
    cmd = "%s sysname=%s -i facts.inp" % (os.path.join(ppdg.CHARMM, 'charmm'), fname)
    out, err, ret = ppdg.tools.execute(cmd)
    with open('%s-facts.out' % (fname), 'w') as fp:
        fp.write(out)
    with open('%s-facts.err' % (fname), 'w') as fp:
        fp.write(err)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))
    for line in out.splitlines():
        if line.startswith('ENER EXTERN>'):
            _, _, vdw, elec, _, _, _ = line.split()
        if line.startswith('ENER FCTPOL>'):
            _, _, gb = line.split()
        if line.startswith('ENER FCTNPL>'):
            _, _, asp = line.split()
    return float(elec), float(vdw), float(gb), float(asp)

def facts(wrkdir):
    """
        Get FACTS solvation term for protein-protein as implemented in CHARMM
    """
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting FACTS scoring...")
    ppdg.link_data('facts.inp')

    cpx_elec, cpx_vdw, cpx_gb, cpx_asp = run_facts('complex-chm')
    rec_elec, rec_vdw, rec_gb, rec_asp = run_facts('receptor-chm')
    lig_elec, lig_vdw, lig_gb, lig_asp = run_facts('ligand-chm')

    desc = dict()
    desc['FACTS_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['FACTS_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    desc['FACTS_GB']   = cpx_gb - (rec_gb + lig_gb)
    desc['FACTS_ASP']  = cpx_asp - (rec_asp + lig_asp)
    time_end = timer()
    desc['>TIME_FACTS'] = time_end - time_start
    os.chdir(basepath)

    return desc


