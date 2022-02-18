#!/usr/bin/env python3

import os
from .modeller import modeller_veryfast, modeller_fast, modeller_slow
from .rosetta import rosetta1
#from .charmify import charmm_model
import ppdx
import logging
log = logging.getLogger(__name__)

def make_model(wrkdir, protocol, tpl_complex, seq_complex): #, tpl_receptor=None, seq_receptor=None, tpl_ligand=None, seq_ligand=None):
    """
        Given sequence, template, working directory and modeling protocol,
        create a 3D atomistic model of the complex.
        Save the output as model.pdb in wrkdir.
    """
    # If the model is already there, just return.
    if os.path.isdir(wrkdir):
        if os.path.isfile(os.path.join(wrkdir, 'model.pdb')):
            return dict()
    # Else, make the model!
    if protocol == 'modeller_veryfast':
        time = modeller_veryfast(seq_complex, tpl_complex, wrkdir)
    elif protocol == 'modeller_fast':
        time = modeller_fast(seq_complex, tpl_complex, wrkdir)
    elif protocol == 'modeller_slow':
        time = modeller_slow(seq_complex, tpl_complex, wrkdir)
    elif protocol == 'rosetta':
        time = rosetta1(seq_complex, tpl_complex, wrkdir)
    else:
        raise ValueError('Unknown requested protocol %s. Valid values are modeller_veryfast, modeller_fast, modeller_slow, and rosetta.' % protocol)
    return {'>TIME_makemodel':time}

def split_complex(wrkdir, nchains):
    """
        Split the model-chm.pdb/cor/psf in wrkdir in a ligand-chm.pdb/cor/psf, 
        receptor-chm.pdb/cor/psf, complex-chm.pdb/cor/psf.
        nchains is a tuple containing the number of chains in the receptor 
        and the number of chains in the ligand.
        Assume the chains are called alphabetically (as done by Modeller).
    """
    basepath = os.getcwd()
    os.chdir(wrkdir)

    if not os.path.isfile('ligand-chm.psf') or not os.path.isfile('receptor-chm.psf') or not os.path.isfile('complex-chm.psf') \
        or not os.path.isfile('ligand.pdb') or not os.path.isfile('receptor.pdb') or not os.path.isfile('complex.pdb') \
        or not os.path.isfile('ligandB.pdb') or not os.path.isfile('receptorA.pdb') or not os.path.isfile('complexAB.pdb'):

        log.info('Splitting model-chm.psf in ligand, receptor and complex.')
        cpx = ppdx.Pdb(os.path.join(wrkdir, 'model-chm.pdb'))
        cpx.segid2chain()
        nchains_tot = len(cpx.split_by_chain())
        lrec = nchains[0]
        llig = nchains[1]
        if llig+lrec!=nchains_tot:
            raise ValueError('PDB %s contains %d chains, but ligand (%d) and receptor (%d) contain %d' % (wrkdir, nchains_tot, llig, lrec, llig+lrec))
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        crec = alphabet[:lrec]
        clig = alphabet[lrec:nchains_tot]
        sele_rec = ' .or. '.join(['segid %s' % (c) for c in crec])
        sele_lig = ' .or. '.join(['segid %s' % (c) for c in clig])
        ppdx.link_data('extract.inp')

        cmd = '%s basename=%s sel="%s" outname=%s ffpath=%s -i extract.inp >%s 2>&1' % (ppdx.CHARMM, 'model-chm', sele_rec, 'receptor-chm', ppdx.FFPATH, 'receptor-chm.out')
        ret = ppdx.tools.execute(cmd)
        if ret!=0:
            raise ValueError("Charmm failed.")

        cmd = '%s basename=%s sel="%s" outname=%s ffpath=%s -i extract.inp >%s 2>&1' % (ppdx.CHARMM, 'model-chm', sele_lig, 'ligand-chm', ppdx.FFPATH, 'ligand-chm.out')
        ret = ppdx.tools.execute(cmd)
        if ret!=0:
            raise ValueError("Charmm failed.")

        for f in ['cor', 'pdb', 'psf']:
            if not os.path.isfile('complex-chm.'+f):
                os.symlink('model-chm.'+f, 'complex-chm.'+f)

        rec = ppdx.Pdb('receptor-chm.pdb')
        lig = ppdx.Pdb('ligand-chm.pdb')
        rec.make_standard()
        lig.make_standard()
        cpx = rec+lig
        rec.write('receptor.pdb')
        lig.write('ligand.pdb')
        cpx.write('complex.pdb')
        rec.set_chain('A')
        rec.set_segid('A')
        lig.set_chain('B')
        lig.set_segid('B')
        cpx = rec+lig
        rec.write('receptorA.pdb')
        lig.write('ligandB.pdb')
        cpx.write('complexAB.pdb')

    os.chdir(basepath)

