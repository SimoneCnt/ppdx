#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

import logging
log = logging.getLogger(__name__)
# Very verbose output
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')
# Only errors
#logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

import ppdx
ppdx.config.cread('config-ppdx.ini')

def main():

    # Set the working directory
    ppdx.WRKDIR= os.path.join(os.getcwd(), "models")

    # Name of the pkl file containing the details of the scoring function
    PKLFILE = 'regression-hiv-rfha-veryfast-6224.pkl'

    # Read sequences
    sequences = ppdx.tools.read_multi_fasta('sequences.seq')
    sequences.update(ppdx.tools.read_multi_fasta('seaman.seq'))

    # Get the list of antigens in the Seaman panel
    SEAMAN = ppdx.tools.read_multi_fasta('seaman.seq').keys()

    # Compute the binding affinity for the VRC01 and VRC01GL antibodies against the BG505 SOSIP antigen
    inputs = list()
    inputs.append(['BG505_VRC01'  , '%s/%s' % (sequences['BG505'], sequences['VRC01'])  , (1,2), 'tpl_5fyj.pdb'])
    inputs.append(['BG505_VRC01GL', '%s/%s' % (sequences['BG505'], sequences['VRC01GL']), (1,2), 'tpl_5fyj.pdb'])
    res = ppdx.eval_pkl(PKLFILE, inputs, nmodels=12)
    print("Binding affinity VRC01       : %8.3f" % (res['BG505_VRC01']))
    print("Binding affinity VRC01GL     : %8.3f" % (res['BG505_VRC01GL']))

    # Make a list of all dg to compute to evaluate the breadth of VRC01 and VRC01GL
    inputs = list()
    for ag in SEAMAN:
        inputs.append(['%s_VRC01'%(ag)  , '%s/%s' % (sequences[ag], sequences['VRC01'])  , (1,2), 'tpl_5fyj.pdb'])
        inputs.append(['%s_VRC01GL'%(ag), '%s/%s' % (sequences[ag], sequences['VRC01GL']), (1,2), 'tpl_5fyj.pdb'])

    # Compute everything
    res = ppdx.eval_pkl(PKLFILE, inputs, nmodels=12, config='pool:12')

    # Recall all computed values and save them into two arrays
    vrc01   = np.array([ res['%s_VRC01'%(ag)] for ag in SEAMAN ])
    vrc01gl = np.array([ res['%s_VRC01GL'%(ag)] for ag in SEAMAN ])

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

