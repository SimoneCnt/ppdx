#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def run_gb(fname, mode):
    if mode not in ['gbsw', 'gbmv']:
        raise ValueError('Mode must be either gbsw or gbmv')
    cmd = "%s sysname=%s -i %s.inp" % (os.path.join(ppdg.CHARMM, 'charmm'), fname, mode)
    out, err, ret = ppdg.tools.execute(cmd)
    with open('%s-%s.out' % (fname, mode), 'w') as fp:
        fp.write(out)
    with open('%s-%s.err' % (fname, mode), 'w') as fp:
        fp.write(err)
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("Charmm failed while running < %s >." % (cmd))
    for line in out.splitlines():
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

    cpx_elec, cpx_vdw, cpx_gb, cpx_asp = run_gb('complex-chm', 'gbsw')
    rec_elec, rec_vdw, rec_gb, rec_asp = run_gb('receptor-chm', 'gbsw')
    lig_elec, lig_vdw, lig_gb, lig_asp = run_gb('ligand-chm', 'gbsw')

    desc = dict()
    desc['GBSW_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc['GBSW_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    desc['GBSW_GB']   = cpx_gb - (rec_gb + lig_gb)
    desc['GBSW_ASP']  = cpx_asp - (rec_asp + lig_asp)
    time_end = timer()
    desc['>TIME_GBSW'] = time_end - time_start
    os.chdir(basepath)

    return desc


