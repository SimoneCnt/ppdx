#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Agg')

import mdtraj

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

import ppdx
ppdx.config.cread('config-ppdx.ini')
ppdx.WRKDIR = os.path.join(os.getcwd(), "models")

def decompress(base_wrkdir, protocol, nmodel):
    wrkdir = os.path.join(base_wrkdir, protocol+'_%s' % (str(nmodel)))
    if os.path.isfile(wrkdir+'.tgz'):
        if os.path.isdir(wrkdir):
            raise ValueError('Both tgz and dir exists for %s' % (wrkdir))
        log.info('Un-compressing %s.tgz' % (wrkdir))
        basepath = os.getcwd()
        os.chdir(base_wrkdir)
        ppdx.tools.execute('tar -xf %s_%s.tgz && rm %s_%s.tgz' % (protocol, str(nmodel), protocol, str(nmodel)))
        os.chdir(basepath)

def compute_rmsd(wrkdir, protocol, nmodels):
    structs = list()
    for i in range(nmodels):
        decompress(wrkdir, protocol, i)
        s = mdtraj.load('%s/%s_%d/model.pdb' % (wrkdir, protocol, i))
        ha = s.topology.select('mass > 4')
        s = s.atom_slice(ha)
        structs.append(s)
    rmsd = list()
    for i in range(0, nmodels):
        for j in range(i+1, nmodels):
            rmsd.append(10.0*mdtraj.rmsd(structs[i], structs[j])[0])
    rmsd = np.array(rmsd)
    rmin = np.amin(rmsd)
    rmax = np.amax(rmsd)
    ravg = np.average(rmsd)
    return rmin, rmax, ravg

def main(protocol):
    print(protocol)
    nmodels = 20
    sequences = ppdx.tools.read_multi_fasta('ppdb/ppdb.seq')
    db1 = list()
    with open('ppdb/ppdb.txt') as fp:
        for line in fp:
            if line[0]=='#':
                continue
            db1.append(line.split())
    # Prepare list of what to compute
    fp = open('rmsd-%s.txt' % (protocol), 'w')
    fp.write('#Name                        rmsd_min  rmsd_max  rmsd_avg\n')
    for name, recn, lign, dgexp, tpl in db1:
        print(name)
        base_wrkdir = os.path.join(ppdx.WRKDIR, name)
        rmin, rmax, ravg = compute_rmsd(base_wrkdir, protocol, nmodels)
        fp.write('%-27s %8.3f  %8.3f  %8.3f\n' % (name, rmin, rmax, ravg))
        ppdx.clean()

main('modeller_fast')
main('rosetta')

