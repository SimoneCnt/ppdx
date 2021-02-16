#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def run_facts(fname):
    outfile = '%s-facts.out' % (fname)
    cmd = "%s sysname=%s domini=0 -i facts.inp >%s 2>&1" % (os.path.join(ppdg.CHARMM, 'charmm'), fname, outfile)
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))
    with open(outfile) as fp:
        for line in fp.readlines():
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

    # Minimize complex
    cmd = "%s sysname=complex-chm domini=1 -i facts.inp >facts_mini.out 2>&1" % (os.path.join(ppdg.CHARMM, 'charmm'))
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))

    # Split minimized complex

    with open('nchains.dat') as fp:
        lrec, llib = (int(v) for v in fp.readline().split())
    nchains_tot = lrec+llib
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    crec = alphabet[:lrec]
    clig = alphabet[lrec:nchains_tot]
    sele_rec = ' .or. '.join(['segid %s' % (c) for c in crec])
    sele_lig = ' .or. '.join(['segid %s' % (c) for c in clig])
    ppdg.link_data('extract.inp')

    cmd = '%s basename=%s sel="%s" outname=%s -i extract.inp >%s 2>&1' % \
            (os.path.join(ppdg.CHARMM, 'charmm'), 'complex-chm-facts', sele_rec, 'receptor-chm-facts', 'receptor-chm-facts.out')
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed.")

    cmd = '%s basename=%s sel="%s" outname=%s -i extract.inp >%s 2>&1' % \
            (os.path.join(ppdg.CHARMM, 'charmm'), 'complex-chm-facts', sele_lig, 'ligand-chm-facts', 'ligand-chm-facts.out')
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed.")


    cpx_elec, cpx_vdw, cpx_gb, cpx_asp = run_facts('complex-chm-facts')
    rec_elec, rec_vdw, rec_gb, rec_asp = run_facts('receptor-chm-facts')
    lig_elec, lig_vdw, lig_gb, lig_asp = run_facts('ligand-chm-facts')

    desc = dict()
    desc['FACTS_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['FACTS_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    desc['FACTS_GB']   = cpx_gb - (rec_gb + lig_gb)
    desc['FACTS_ASP']  = cpx_asp - (rec_asp + lig_asp)
    time_end = timer()
    desc['>TIME_FACTS'] = time_end - time_start
    os.chdir(basepath)

    return desc


