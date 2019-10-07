#!/usr/bin/env python3

import os
import ppdg
import logging
log = logging.getLogger(__name__)

def solvate(fname, align=False, margin=9, conc=0.15):
    '''
        Suppose fname.pdb and fname.psf exist!
    '''

    basepath = os.getcwd()
    wrkdir, name = os.path.split(fname)

    if os.path.isfile(os.path.join(wrkdir, '%s-ion.psf' % (name))):
        log.info('Already present!')
        return
    else:
        log.info("Solvating %s" % (name))

    if align:
        nchains = len(ppdg.Pdb('%s.pdb' % (name)).split_by_chain())
        if nchains==1:
            align=False
        elif nchains==2:
            rec = 'A'
            lig = 'B'
        elif nchains==3:
            rec = 'A'
            lig = 'B C'
        elif nchains==4:
            rec = 'A B'
            lig = 'C D'
        else:
            raise ValueError('Unexpected number of chains %d with align. (Max supported: 4)' % (nchains))

    if len(wrkdir)>0:
        os.chdir(wrkdir)

    with open('vmd.tcl', 'w') as fp:
        fp.write('mol new %s.psf\n' % (name))
        if os.path.isfile('%s.pdb' % (name)):
            fp.write('mol addfile %s.pdb\n' % (name))
        elif os.path.isfile('%s.cor' % (name)):
            fp.write('mol addfile %s.cor\n' % (name))
        else:
            raise ValueError('Impossible to find neither %s.pdb nor %s.cor' % (fname, fname))
        fp.write('set all [atomselect top "all"]\n')
        fp.write('set com [measure center $all]\n')
        fp.write('$all moveby [vecscale -1.0 $com]\n')
        if align:
            fp.write('set lig [atomselect top "segid %s and name CA"]\n' % (lig))
            fp.write('set rec [atomselect top "segid %s and name CA"]\n' % (rec))
            fp.write('set com1 [measure center $lig]\n')
            fp.write('set com2 [measure center $rec]\n')
            fp.write('set matrix [transvecinv $com2]\n')
            fp.write('$all move $matrix\n')
        fp.write('$all writepdb "%s-align.pdb"\n' % (name))
        fp.write('package require solvate\n')
        fp.write('solvate %s.psf %s-align.pdb -t %f -o %s-solv\n' % (name, name, margin, name))
        fp.write('package require autoionize\n')
        fp.write('autoionize -psf %s-solv.psf -pdb %s-solv.pdb -sc %f -cation POT -o %s-ion\n' % (name, name, conc, name))
        fp.write('quit\n')

    ret = ppdg.tools.execute('vmd -dispdev text -e vmd.tcl')

    os.chdir(basepath)
    if ret!=0:
        raise ValueError("VMD failed.")


