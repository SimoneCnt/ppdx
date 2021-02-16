#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

import logging
log = logging.getLogger(__name__)
# Very verbose output
#logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')
# Only errors
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

import ppdg
ppdg.config.cread('config-ppdg.ini')

def compute(abname, agname, pkl, template='tpl_5fyj.pdb', nmodels=12, ncores=1):
    """
        Utility function to prepare the inputs for ppdg.eval_pkl
        pkl: contains the scoring function
        template: a pdb structure to use as a template
        nmodels: the number of models to create and comput averages (usually between 5 and 15 are good numbers)
        ncores: number of cores to run in parallel
    """
    # Gather the sequences of the antibody and the antigen. We need a single sequence with "/" as chain separator.
    abname = abname.upper()
    agname = agname.upper()
    sequences = ppdg.tools.read_multi_fasta('sequences.seq')
    sequences.update(ppdg.tools.read_multi_fasta('seaman.seq'))
    abseq = sequences[abname]
    agseq = sequences[agname]
    cpxseq = '%s/%s' % (agseq, abseq)
    # Number of chains in the receptor (the antigen) and ligand (the antibody)
    nchains = (1, 2)
    # A custom name, whatever you like
    name = '%s__%s' % (agname, abname)
    name = name.lower()
    # Working directory for this protein-protein complex
    wrkdir = os.path.join(ppdg.WRKDIR, name)
    # Compute the binding affinity
    # Specify ncores=n if you want to run in parallel using n cores.
    dg = ppdg.eval_pkl(pkl, cpxseq, nchains, template, name=name, wrkdir=wrkdir, nmodels=nmodels, ncores=ncores)
    return dg

def main():

    # Set the working directory
    ppdg.WRKDIR= os.path.join(os.getcwd(), "models")

    # Name of the pkl file containing the details of the scoring function
    PKLFILE = 'regression-hiv-rfha-veryfast-6224.pkl'

    # Get the list of antigens in the Seaman panel
    SEAMAN = ppdg.tools.read_multi_fasta('seaman.seq').keys()

    # Compute the binding affinity for the VRC01 and VRC01GL antibodies against the BG505 SOSIP antigen
    dg_vrc01   = compute('VRC01',   'BG505', PKLFILE, 'tpl_5fyj.pdb', ncores=6)
    dg_vrc01gl = compute('VRC01GL', 'BG505', PKLFILE, 'tpl_5fyj.pdb', ncores=6)
    print("Binding affinity VRC01       : %8.3f" % (dg_vrc01))
    print("Binding affinity VRC01GL     : %8.3f" % (dg_vrc01gl))

    # Make a list of all dg to compute to evaluate the breadth of VRC01 and VRC01GL
    to_compute = list()
    for ag in SEAMAN:
        to_compute.append(['VRC01',   ag, PKLFILE, 'tpl_5fyj.pdb'])
        to_compute.append(['VRC01GL', ag, PKLFILE, 'tpl_5fyj.pdb'])

    # Compute everything in parallel (using 10 cores here)
    with Pool(10) as p:
        p.starmap(compute, to_compute)

    # Recall all computed values and save them into two arrays
    vrc01   = np.array([ compute('VRC01',   ag, PKLFILE, 'tpl_5fyj.pdb') for ag in SEAMAN ]).ravel()
    vrc01gl = np.array([ compute('VRC01GL', ag, PKLFILE, 'tpl_5fyj.pdb') for ag in SEAMAN ]).ravel()

    # Compute the breadth. Use a cutoff energy of -9.6
    b_vrc01   = np.sum(vrc01<-9.6)/len(SEAMAN)
    b_vrc01gl = np.sum(vrc01gl<-9.6)/len(SEAMAN)
    print("Computed breadth for VRC01   : %8.3f" % (b_vrc01))
    print("Computed breadth for VRC01GL : %8.3f" % (b_vrc01gl))

    # Plot the histogram of the binding affinities for the two antibodies against the antigens in the Seaman panel
    plt.style.use('seaborn-deep')
    plt.hist([vrc01, vrc01gl], bins='auto', label=['VRC01 (B=%.2f)'%(b_vrc01),'VRC01GL (B=%.2f)'%(b_vrc01gl)], edgecolor='k')
    plt.legend()
    plt.xlabel('Binding Affinity')
    plt.ylabel('Count')
    plt.savefig('histogram_breadth.png', bbox_inches='tight', dpi=300)
    #plt.show()

main()

