#!/usr/bin/env python3

import os
import json

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

import ppdg
ppdg.config.cread('config-ppdg.ini')


def save_descriptors_json(inputs, jname=None):
    def getone(n):
        wrkdir = os.path.join(ppdg.WRKDIR, n)
        descfile = os.path.join(wrkdir, 'descriptors.json')
        if not os.path.isfile(descfile):
            raise ValueError("Missing file descriptors.json in directory %s\nBe sure to have called get_descriptors before!" % (wrkdir))
        with open(descfile, 'r') as fp:
            desc = { n : json.load(fp) }
        return desc

    if type(inputs) is str:
        res = getone(inputs)
    elif type(inputs) is list:
        res = dict()
        for name, _, _, _ in inputs:
            res.update(getone(name))
    else:
        raise ValueError("Cound not understand inputs. Should be a string or a list.")
    if jname:
        with open(jname, 'w') as fp:
            json.dump(res, fp, indent=4, sort_keys=True)
    return res


def compute(dbpath, nmodels=12, config='pool:12', protocol='modeller_fast'):

    desclist = list()
    desclist += ['HB_BH', 'HB_WN', 'HB_KS']
    desclist += ['BSA', 'BSA_C', 'BSA_A', 'BSA_P', 'NIS_P', 'NIS_C', 'NIS_A', 'NRES']
    desclist += ['sticky_tot', 'sticky_avg']
    desclist += ['IC_TOT', 'IC_AA', 'IC_PP', 'IC_CC', 'IC_AP', 'IC_CP', 'IC_AC']
    desclist += ['ZRANK', 'ZRANK2']
    desclist += ['pyDock', 'pyDock_elec', 'pyDock_vdw', 'pyDock_desolv']
    #desclist += ['ATTRACT']
    #desclist += ['FireDock', 'FireDock_aVdW', 'FireDock_rVdW', 'FireDock_ACE', 'FireDock_inside',
    #            'FireDock_aElec', 'FireDock_rElec', 'FireDock_laElec', 'FireDock_lrElec',
    #            'FireDock_hb', 'FireDock_piS', 'FireDock_catpiS', 'FireDock_aliph']
    desclist += ['RF_HA_SRS', 'RF_CB_SRS_OD']
    desclist += ['ipot_aace167', 'ipot_aace18', 'ipot_aace20', 'ipot_rrce20']
    #desclist += ['SOAP-PP-Pair', 'SOAP-Protein-OD']
    desclist += ['DOPE', 'DOPE-HR']
    desclist += ['AGBNP']
    desclist += ['FACTS_ELEC', 'FACTS_VDW', 'FACTS_GB', 'FACTS_ASP']
    desclist += ['GBMV_ELEC', 'GBMV_VDW', 'GBMV_GB', 'GBMV_ASP']
    desclist += ['GBSW_ELEC', 'GBSW_VDW', 'GBSW_GB', 'GBSW_ASP']
    desclist += ['CDIE_ELEC', 'CDIE_VDW']
    desclist += ['RDIE_ELEC', 'RDIE_VDW']
    #desclist += ['OMM_vacuum', 'OMM_HCT', 'OMM_OBC1', 'OMM_OBC2', 'OMM_GBn', 'OMM_GBn2']
    desclist += ['FoldX', 'FoldX_backbone_hbond', 'FoldX_sidechain_hbond', 'FoldX_vdw', 
                'FoldX_elec', 'FoldX_solvation_polar', 'FoldX_solvation_hydrophobic', 
                'FoldX_entropy_sidechain', 'FoldX_entropy_mainchain']
    desclist += ['Rosetta_dg']
    desclist += ['Prodigy_IC_NIS']

    desclist += ['ENM_EXP']
    desclist += ['ENM_R6']


    # Prepare list of what to compute
    inputs = list()
    sequences = ppdg.tools.read_multi_fasta('ppdb/ppdb.seq')
    with open('ppdb/ppdb.txt') as fp:
        for line in fp:
            if line[0]=='#':
                continue
            name, recn, lign, dgexp, tpl = line.split()
            sequence = sequences[name]
            template = '%s/%s' % (dbpath, tpl)
            nchains = (int(recn), int(lign))
            inputs.append([name, sequence, nchains, template])

    # Compute the descriptors
    ppdg.eval_descriptors(protocol, desclist, inputs, nmodels=nmodels, config=config)
    save_descriptors_json(inputs, 'descriptors-all.json')


if __name__=='__main__':
    ppdg.WRKDIR = os.path.join(os.getcwd(), "models")
    for n in range(20):
        for protocol in ['modeller_fast', 'modeller_veryfast', 'modeller_slow' ]: #, 'rosetta']:
            print(n, protocol)
            compute('ppdb/', nmodels=n+1, config='pool:12', protocol=protocol)
            ppdg.clean()
            if os.path.isfile('kill'):
                print('Kill!')
                quit()


