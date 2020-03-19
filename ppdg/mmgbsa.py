#!/usr/bin/env python3

import os, json
import hashlib
import joblib
import numpy as np
import simtk.unit as unit
from multiprocessing import Pool
import ppdg
import logging
log = logging.getLogger(__name__)

def mmgbsa(wrkdir, template, sequence, nchains):
    """
        Main function!
        Given:
            - The working directory to use
            - A PDB file to use as template for the modelling
            - The sequence of the whole protein-protein complex
            - A tuple indicating how many chains are in the receptor and in the ligand [ nchains = ( nchains_rec, nchains_lig) ]
        This function
            - Makes the model
            - Runs MD
            - Computes mmgbsa
    """

    # Check template file exists
    if not os.path.isfile(template):
        raise ValueError('Template file %s does not exist!' % (template))

    # Check nchains. A tuple. First is the number of chains in the receptor, 
    # second the number of chains in the ligand.
    if len(nchains)!=2:
        raise ValueError('Cannot understand nchains = %s. Use (nchains_receptor, nchains_ligand).' % (nchains))
    nchains_rec = nchains[0]
    nchains_lig = nchains[1]
    nchains_tot = len(sequence.split('/'))
    if nchains_rec + nchains_lig != nchains_tot:
        raise ValueError('You said the receptor has %d chains and the ligand %d, but the sequence has a total of %d chains.' 
                % (nchains_rec, nchains_lig, nchains_tot))

    # Create model
    scores = ppdg.makemodel.make_model(wrkdir, 'modeller_veryfast', template, sequence)
    ppdg.makemodel.charmify(os.path.join(wrkdir, 'model.pdb'), nsteps=10)
    ppdg.makemodel.split_complex(wrkdir, nchains)

    # Solvate
    basepath = os.getcwd()
    os.chdir(wrkdir)
    ppdg.makemodel.solvate_chm('complex-chm', buf=10.0)

    # Run MD complex
    for i in range(1):
        outname = 'md_complex_%d' % (i)
        if i>0:
            prevxml = 'md_complex_%d-final.xml' % (i-1)
        else:
            prevxml = None
        if os.path.isfile('md_complex_%d-final.xml' % (i)):
            continue
        log.info('Running MD for complex for step %d/%d (outname=%s, prevxml=%s)' % (i+1, 3, outname, prevxml))
        ppdg.runomm(
            psffile='complex-chm_center.psf',
            crdfile='complex-chm_center.cor',
            simultime=10*unit.nanoseconds,
            cutoff=10*unit.angstrom,
            outname=outname,
            prevxml=prevxml
        )

    os.chdir(basepath)


