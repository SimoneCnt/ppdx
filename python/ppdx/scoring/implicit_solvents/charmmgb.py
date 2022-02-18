#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import ppdx
import logging
log = logging.getLogger(__name__)

def run_gb(fname, mode):
    outfile = '%s-%s.out' % (fname, mode)
    cmd = "%s sysname=%s domini=0 ffpath=%s -i %s.inp >%s 2>&1" % (ppdx.CHARMM, fname, ppdx.FFPATH, mode, outfile)
    ret = ppdx.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s >." % (cmd))
    if mode in ['gbsw', 'gbmv']:
        with open(outfile) as fp:
            for line in fp.readlines():
                if line.startswith('ENER EXTERN>'):
                    _, _, vdw, elec, _, asp, _ = line.split()
                if line.startswith('ENER PBEQ>'):
                    _, _, _, _, gb, _, _ = line.split()
    elif mode in ['facts']:
        with open(outfile) as fp:
            for line in fp.readlines():
                if line.startswith('ENER EXTERN>'):
                    _, _, vdw, elec, _, _, _ = line.split()
                if line.startswith('ENER FCTPOL>'):
                    _, _, gb = line.split()
                if line.startswith('ENER FCTNPL>'):
                    _, _, asp = line.split()
    elif mode in ['cdie', 'rdie']:
        gb = 0.0
        asp = 0.0
        with open(outfile) as fp:
            for line in fp.readlines():
                if line.startswith('ENER EXTERN>'):
                    _, _, vdw, elec, _, _, _ = line.split()
    else:
        raise ValueError('Cannot understand mode %s' % (mode))

    return float(elec), float(vdw), float(gb), float(asp)


def charmm_gen(wrkdir, lmode):
    """
        Get CHARMM solvation term for protein-protein.
    """
    umode = lmode.upper()
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting %s scoring..." % (umode))
    ppdx.link_data(lmode+'.inp')

    # Minimize complex
    cmd = "%s sysname=complex-chm domini=1 ffpath=%s -i %s.inp >%s_mini.out 2>&1" % (ppdx.CHARMM, ppdx.FFPATH, lmode, lmode)
    ret = ppdx.tools.execute(cmd)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s > in %s" % (cmd, wrkdir))
    os.remove('complex-chm-%s.psf' % (lmode))
    os.symlink('complex-chm.psf', 'complex-chm-%s.psf' % (lmode))

    # Split minimized complex
    with open('nchains.dat') as fp:
        lrec, llib = (int(v) for v in fp.readline().split())
    nchains_tot = lrec+llib
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    crec = alphabet[:lrec]
    clig = alphabet[lrec:nchains_tot]
    sele = dict(
            receptor = ' .or. '.join(['segid %s' % (c) for c in crec]),
            ligand = ' .or. '.join(['segid %s' % (c) for c in clig])
        )
    ppdx.link_data('extract.inp')

    # Minimize the receptor and the ligand
    for sysname in ['receptor', 'ligand']:
        cmd = '%s basename=%s sel="%s" outname=%s ffpath=%s -i extract.inp >%s 2>&1' % \
            (ppdx.CHARMM, 'complex-chm-'+lmode, sele[sysname], sysname+'-chm-'+lmode, ppdx.FFPATH, sysname+'-chm-'+lmode+'.out')
        ret = ppdx.tools.execute(cmd)
        if ret!=0:
            raise ValueError("Charmm failed while running < %s > in %s" % (cmd, wrkdir))
        os.remove('%s-chm-%s.psf' % (sysname, lmode))
        os.symlink('%s-chm.psf' % (sysname), '%s-chm-%s.psf' % (sysname, lmode))

    # Get the energies
    cpx_elec, cpx_vdw, cpx_gb, cpx_asp = run_gb('complex-chm', lmode)
    rec_elec, rec_vdw, rec_gb, rec_asp = run_gb('receptor-chm', lmode)
    lig_elec, lig_vdw, lig_gb, lig_asp = run_gb('ligand-chm', lmode)

    desc = dict()
    desc[umode+'_ELEC'] = cpx_elec - (rec_elec + lig_elec)
    desc[umode+'_VDW']  = cpx_vdw - (rec_vdw + lig_vdw)
    if lmode in ['facts', 'gbmv', 'gbsw']:
        desc[umode+'_GB']   = cpx_gb - (rec_gb + lig_gb)
        desc[umode+'_ASP']  = cpx_asp - (rec_asp + lig_asp)
        desc[umode+'_POL'] = desc[umode+'_ELEC'] + desc[umode+'_GB']
        desc[umode+'_TOT'] = desc[umode+'_POL'] + desc[umode+'_VDW']
    else:
        desc[umode+'_TOT'] = desc[umode+'_ELEC'] + desc[umode+'_VDW']
    time_end = timer()
    desc['>TIME_'+umode] = time_end - time_start
    os.chdir(basepath)

    return desc


def gbmv(wrkdir):
    return charmm_gen(wrkdir, 'gbmv')

def gbsw(wrkdir):
    return charmm_gen(wrkdir, 'gbsw')

def facts(wrkdir):
    return charmm_gen(wrkdir, 'facts')

def cdie(wrkdir):
    return charmm_gen(wrkdir, 'cdie')

def rdie(wrkdir):
    return charmm_gen(wrkdir, 'rdie')

