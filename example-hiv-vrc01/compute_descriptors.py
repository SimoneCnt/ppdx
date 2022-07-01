#!/usr/bin/env python3

import os
import json

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

import ppdx
ppdx.config.cread('config-ppdx.ini')

def compute(dbpath, nmodels=12, config='pool:12', protocol='modeller_fast'):

    # Create list of descriptors to compute
    desc = list()
    # Docking
    desc += ['ZRANK', 'ZRANK2']
    desc += ['pyDock']
    desc += ['ATTRACT']
    desc += ['FireDock']
    # Statistical Potentials
    desc += ['RF_HA_SRS', 'RF_CB_SRS_OD']
    desc += ['ipot_aace167', 'ipot_aace18', 'ipot_aace20', 'ipot_rrce20']
    desc += ['DOPE', 'DOPE-HR']
    # Implicit Solvents
    desc += ['FACTS_ELEC', 'FACTS_VDW', 'FACTS_GB', 'FACTS_ASP', 'FACTS_POL', 'FACTS_TOT']
    desc += ['GBSW_ELEC', 'GBSW_VDW', 'GBSW_GB', 'GBSW_ASP', 'GBSW_POL', 'GBSW_TOT']
    # Folding
    desc += ['FoldX']
    desc += ['Rosetta_dg']
    # Binding
    desc += ['Prodigy_IC_NIS']

    # Prepare list of what to compute
    inputs = list()
    sequences = ppdx.tools.read_multi_fasta('ppdb/ppdb.seq')
    with open('ppdb/ppdb.txt') as fp:
        for line in fp:
            if line[0]=='#':
                continue
            name, recn, lign, dgexp, tpl = line.split()
            sequence = sequences[name.upper()]
            template = os.path.join(dbpath, tpl)
            nchains = (int(recn), int(lign))
            inputs.append([name, sequence, nchains, template])

    # Compute the descriptors
    ppdx.eval_descriptors(protocol, desc, inputs, nmodels=nmodels, config=config)
    ppdx.save_descriptors_json(inputs, 'descriptors-all.json')


if __name__=='__main__':
    ppdx.WRKDIR = os.path.join(os.getcwd(), "models")
    for n in range(20):
        for protocol in ['modeller_veryfast', 'modeller_fast', 'modeller_slow', 'rosetta']:
            compute(os.path.join(os.getcwd(), 'ppdb'), nmodels=n+1, config='pool', protocol=protocol)
            ppdx.clean()
            if os.path.isfile('kill'):
                print('Kill!')
                quit()


