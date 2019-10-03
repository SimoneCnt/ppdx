#!/usr/bin/env python3

import os
from .modeller import modeller_veryfast, modeller_fast, modeller_slow
#from .charmify import charmm_model
import ppdg
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
    #if protocol == 'charmm':
    #    time = charmm_model(wrkdir, tpl_complex, tpl_receptor, tpl_ligand)
    #    return time
    else:
        raise ValueError('Unknown requested protocol %s. Valid values are modeller_fast and modeller_slow.' % protocol)

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
        cpx = ppdg.Pdb(os.path.join(wrkdir, 'model-chm.pdb'))
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
        ppdg.link_data('extract.inp')

        cmd = '%s basename=%s sel="%s" outname=%s -i extract.inp' % (os.path.join(ppdg.CHARMM, 'charmm'), 'model-chm', sele_rec, 'receptor-chm')
        out, err, ret = ppdg.tools.execute(cmd)
        with open('receptor-chm.out', 'w') as fp:
            fp.write(out)
        with open('receptor-chm.err', 'w') as fp:
            fp.write(err)
        if ret!=0:
            raise ValueError("Charmm failed.")

        cmd = '%s basename=%s sel="%s" outname=%s -i extract.inp' % (os.path.join(ppdg.CHARMM, 'charmm'), 'model-chm', sele_lig, 'ligand-chm')
        out, err, ret = ppdg.tools.execute(cmd)
        with open('ligand-chm.out', 'w') as fp:
            fp.write(out)
        with open('ligand-chm.err', 'w') as fp:
            fp.write(err)
        if ret!=0:
            raise ValueError("Charmm failed.")

        for f in ['cor', 'pdb', 'psf']:
            if not os.path.isfile('complex-chm.'+f):
                os.symlink('model-chm.'+f, 'complex-chm.'+f)

        rec = ppdg.Pdb('receptor-chm.pdb')
        lig = ppdg.Pdb('ligand-chm.pdb')
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

        

#def split_model(wrkdir, nchains):
#    """
#        Split the model.pdb in wrkdir in a ligand.pdb, receptor.pdb, 
#        complex.pdb and complexAB.pdb. The last contains the whole receptor 
#        as chain A and the full ligand as chain B. Useful for ZRANK.
#        nchains is a tuple containing the number of chains in the receptor 
#        and the number of chains in the ligand.
#    """
#    basepath = os.getcwd()
#    os.chdir(wrkdir)
#    if not os.path.isfile('ligand.pdb') or not os.path.isfile('receptor.pdb') or not os.path.isfile('complex.pdb') or not os.path.isfile('complexAB.pdb'):
#        log.info('Splitting model.pdb in ligand, receptor and complex.')
#        cpx = ppdg.Pdb(os.path.join(wrkdir, 'model-chm.pdb'))
#        cpx.set_beta(1.0)
#        cpx.set_occupancy(1.0)
#        nchains_tot = len(cpx.split_by_chain())
#        lrec = nchains[0]
#        llig = nchains[1]
#        if llig+lrec!=nchains_tot:
#            raise ValueError('PDB %s contains %d chains, but ligand (%d) and receptor (%d) contain %d' % (wrkdir, nchains_tot, llig, lrec, llig+lrec))
#        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
#        crec = alphabet[:lrec]
#        clig = alphabet[lrec:nchains_tot]
#        rec = cpx.extract(chain=crec)
#        lig = cpx.extract(chain=clig)
#        cpx = rec+lig
#        lig.write('ligand.pdb')
#        rec.write('receptor.pdb')
#        cpx.write('complex.pdb')
#        recA = rec
#        recA.set_chain('A')
#        ligB = lig
#        ligB.set_chain('B')
#        cpxAB = recA + ligB
#        cpxAB.write('complexAB.pdb')
#    os.chdir(basepath)

