#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def run_gb(fname, mode):
    if mode not in ['gbsw', 'gbmv']:
        raise ValueError('Mode must be either gbsw or gbmv')
    outfile = '%s-%s.out' % (fname, mode)
    cmd = "%s sysname=%s domini=0 -i %s.inp >%s 2>&1" % (os.path.join(ppdg.CHARMM, 'charmm'), fname, mode, outfile)
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))
    with open(outfile) as fp:
        for line in fp.readlines():
            if line.startswith('ENER EXTERN>'):
                _, _, vdw, elec, _, asp, _ = line.split()
            if line.startswith('ENER PBEQ>'):
                _, _, _, _, gb, _, _ = line.split()
    return float(elec), float(vdw), float(gb), float(asp)

def gbmv(wrkdir):
    """
        Get GBMV solvation term for protein-protein as implemented in CHARMM
    """
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting GBMV scoring...")
    ppdg.link_data('gbmv.inp')

    # Minimize complex
    cmd = "%s sysname=complex-chm domini=1 -i gbmv.inp >gbmv_mini.out 2>&1" % (os.path.join(ppdg.CHARMM, 'charmm'))
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))

    os.remove('complex-chm-gbmv.psf')
    os.symlink('complex-chm.psf', 'complex-chm-gbmv.psf')

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
            (os.path.join(ppdg.CHARMM, 'charmm'), 'complex-chm-gbmv', sele_rec, 'receptor-chm-gbmv', 'receptor-chm-gbmv.out')
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed.")

    os.remove('receptor-chm-gbmv.psf')
    os.symlink('receptor-chm.psf', 'receptor-chm-gbmv.psf')

    cmd = '%s basename=%s sel="%s" outname=%s -i extract.inp >%s 2>&1' % \
            (os.path.join(ppdg.CHARMM, 'charmm'), 'complex-chm-gbmv', sele_lig, 'ligand-chm-gbmv', 'ligand-chm-gbmv.out')
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed.")

    os.remove('ligand-chm-gbmv.psf')
    os.symlink('ligand-chm.psf', 'ligand-chm-gbmv.psf')

    cpx_elec, cpx_vdw, cpx_gb, cpx_asp = run_gb('complex-chm', 'gbmv')
    rec_elec, rec_vdw, rec_gb, rec_asp = run_gb('receptor-chm', 'gbmv')
    lig_elec, lig_vdw, lig_gb, lig_asp = run_gb('ligand-chm', 'gbmv')

    desc = dict()
    desc['GBMV_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['GBMV_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    desc['GBMV_GB']   = cpx_gb - (rec_gb + lig_gb)
    desc['GBMV_ASP']  = cpx_asp - (rec_asp + lig_asp)
    time_end = timer()
    desc['>TIME_GBMV'] = time_end - time_start
    os.chdir(basepath)

    return desc


def gbsw(wrkdir):
    """
        Get GBSW solvation term for protein-protein as implemented in CHARMM
    """
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting GBSW scoring...")
    ppdg.link_data('gbsw.inp')

    # Minimize complex
    cmd = "%s sysname=complex-chm domini=1 -i gbsw.inp >gbsw_mini.out 2>&1" % (os.path.join(ppdg.CHARMM, 'charmm'))
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
            (os.path.join(ppdg.CHARMM, 'charmm'), 'complex-chm-gbsw', sele_rec, 'receptor-chm-gbsw', 'receptor-chm-gbsw.out')
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed.")

    cmd = '%s basename=%s sel="%s" outname=%s -i extract.inp >%s 2>&1' % \
            (os.path.join(ppdg.CHARMM, 'charmm'), 'complex-chm-gbsw', sele_lig, 'ligand-chm-gbsw', 'ligand-chm-gbsw.out')
    ret = ppdg.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed.")

    cpx_elec, cpx_vdw, cpx_gb, cpx_asp = run_gb('complex-chm-gbsw', 'gbsw')
    rec_elec, rec_vdw, rec_gb, rec_asp = run_gb('receptor-chm-gbsw', 'gbsw')
    lig_elec, lig_vdw, lig_gb, lig_asp = run_gb('ligand-chm-gbsw', 'gbsw')

    desc = dict()
    desc['GBSW_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['GBSW_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    desc['GBSW_GB']   = cpx_gb - (rec_gb + lig_gb)
    desc['GBSW_ASP']  = cpx_asp - (rec_asp + lig_asp)
    time_end = timer()
    desc['>TIME_GBSW'] = time_end - time_start
    os.chdir(basepath)

    return desc


