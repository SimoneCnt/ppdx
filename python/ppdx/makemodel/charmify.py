#!/usr/bin/env python3

import os, shutil
import ppdx
import logging
log = logging.getLogger(__name__)

def charmify(fname, nsteps=100):

    basepath = os.getcwd()
    wrkdir, name = os.path.split(fname)

    basename = ''.join(name.split('.')[0:-1]) + '-chm'

    if os.path.isfile(os.path.join(wrkdir, basename+'.psf')) and os.path.isfile(os.path.join(wrkdir, basename+'.pdb')):
        return
    else:
        log.info("Charmify-ing pdb %s" % (fname))

    os.chdir(wrkdir)

    pdb = ppdx.Pdb(name)
    pdb.fix4charmm()
    pdb.chain2segid()
    pdb.set_occupancy(1.0)
    pdb.set_beta(1.0)
    pdb.remove_hydrogens()

    chains = pdb.split_by_chain()
    nchains = len(chains)
    cmd = ppdx.CHARMM
    cmd += ' nc=%d ' % nchains
    i=1
    for ch, pdb in chains.items():
        pdb.write("chain_%s.pdb" % (ch.lower()))
        cmd += 'c%d=%s ' % (i, ch)
        i += 1
    cmd += 'name=chain_ out=%s ' % (basename)
    cmd += 'nsteps=%d ' % (nsteps)
    cmd += 'ffpath=%s ' % (ppdx.FFPATH)
    cmd += '-i buildgen.inp >%s 2>&1' % (basename+'.out')

    ppdx.link_data('buildgen.inp')
    ppdx.link_data('disu.str')

    ret = ppdx.tools.execute(cmd)
    os.chdir(basepath)
    if ret!=0:
        raise ValueError("Charmm failed while running < %s > in %s" % (cmd, wrkdir))

