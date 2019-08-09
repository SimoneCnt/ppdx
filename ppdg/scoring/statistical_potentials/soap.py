#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import modeller
import modeller.soap_pp
import modeller.soap_protein_od
import ppdg
import logging
log = logging.getLogger(__name__)

def soap_pp(wrkdir):
    """
        Get SOAP-PP Pair score for protein-protein. PDB must contain only two 
        chains. Do not compute the SOAP-PP Atom score.
        [1] G. Q. Dong, H. Fan, D. Schneidman-Duhovny, B. Webb, and A. Sali, 
            "Optimized atomic statistical potentials: assessment of protein 
            interfaces and loops", Bioinformatics, vol. 29, no. 24, 
            pp. 3158-3166, 2013.
    """
    time_start = timer()
    log.info("Getting SOAP-PP-Pair scoring...")
    cpx = os.path.join(wrkdir, 'complexAB.pdb')

    with open(os.path.join(wrkdir, "soap_pp_pair.out"), "w") as fp:
        _stdout = sys.stdout
        sys.stdout = fp
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, cpx)
        cpx = modeller.selection(mdl.chains[:])
        soap_pp_pair = modeller.soap_pp.PairScorer()
        score = cpx.assess(soap_pp_pair)
        sys.stdout.flush()
        sys.stdout = _stdout

    desc = dict()
    desc['SOAP-PP-Pair'] = score
    time_end = timer()
    desc['>TIME_SOAP-PP-Pair'] = time_end - time_start
    return desc


ppdg._soap_scorer = False

def soap_protein_od(wrkdir):
    """
        Get SOAP-Protein score for protein-protein.
        PDB must contain only two chains.
        [1] G. Q. Dong, H. Fan, D. Schneidman-Duhovny, B. Webb, and A. Sali, 
            "Optimized atomic statistical potentials: assessment of protein 
            interfaces and loops", Bioinformatics, vol. 29, no. 24, 
            pp. 3158-3166, 2013.
    """
    if ppdg._soap_scorer==False:
        log.info('Reading SOAP Protein OD Scorer...')
        scorer = modeller.soap_protein_od.Scorer()
        ppdg._soap_scorer = scorer

    time_start = timer()
    log.info("Getting SOAP-Protein-OD scoring...")
    cpx = os.path.join(wrkdir, 'complexAB.pdb')

    with open(os.path.join(wrkdir, "soap_protein_od.out"), "w") as fp:
        _stdout = sys.stdout
        sys.stdout = fp
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, cpx)
        cpx = modeller.selection(mdl.chains[:])
        lig = modeller.selection(mdl.chains[0])
        rec = modeller.selection(mdl.chains[1])
        soap_protein_od = ppdg._soap_scorer
        score = cpx.assess(soap_protein_od) - rec.assess(soap_protein_od) - lig.assess(soap_protein_od)
        sys.stdout.flush()
        sys.stdout = _stdout

    desc = dict()
    desc['SOAP-Protein-OD'] = score
    time_end = timer()
    desc['>TIME_SOAP-Protein-OD'] = time_end - time_start
    return desc
    

