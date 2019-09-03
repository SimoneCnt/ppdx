#!/usr/bin/env python3

import logging
log = logging.getLogger(__name__)

from . import molecular
from . import docking
from . import statistical_potentials
from . import implicit_solvents
from . import folding
from . import entropy
from . import binding

def all_descriptors():
    desc  = list()
    # Molecular
    desc += ['HB_BH', 'HB_WN', 'HB_KS']
    desc += ['BSA', 'BSA_C', 'BSA_A', 'BSA_P', 'NIS_P', 'NIS_C', 'NIS_A', 'NRES']
    desc += ['IC_TOT', 'IC_AA', 'IC_PP', 'IC_CC', 'IC_AP', 'IC_CP', 'IC_AC']
    # Docking
    desc += ['ZRANK', 'ZRANK2']
    desc += ['pyDock', 'pyDock_elec', 'pyDock_vdw', 'pyDock_desolv']
    desc += ['ATTRACT']
    desc += ['FireDock', 'FireDock_aVdW', 'FireDock_rVdW', 'FireDock_ACE', 'FireDock_inside',
                'FireDock_aElec', 'FireDock_rElec', 'FireDock_laElec', 'FireDock_lrElec',
                'FireDock_hb', 'FireDock_piS', 'FireDock_catpiS', 'FireDock_aliph']
    # Statistical Potentials
    desc += ['RF_HA_SRS', 'RF_CB_SRS_OD']
    desc += ['ipot_aace167', 'ipot_aace18', 'ipot_aace20', 'ipot_rrce20']
    desc += ['SOAP-PP-Pair', 'SOAP-Protein-OD']
    desc += ['DOPE', 'DOPE-HR']
    # Implicit Solvents
    desc += ['AGBNP']
    desc += ['FACTS_ELEC', 'FACTS_VDW', 'FACTS_GB', 'FACTS_ASP']
    desc += ['GBMV_ELEC', 'GBMV_VDW', 'GBMV_GB', 'GBMV_ASP']
    desc += ['GBSW_ELEC', 'GBSW_VDW', 'GBSW_GB', 'GBSW_ASP']
    desc += ['CDIE_ELEC', 'CDIE_VDW']
    desc += ['RDIE_ELEC', 'RDIE_VDW']
    desc += ['OMM_vacuum', 'OMM_HCT', 'OMM_OBC1', 'OMM_OBC2', 'OMM_GBn', 'OMM_GBn2']
    # Folding
    desc += ['FoldX', 'FoldX_backbone_hbond', 'FoldX_sidechain_hbond', 'FoldX_vdw', 
                'FoldX_elec', 'FoldX_solvation_polar', 'FoldX_solvation_hydrophobic', 
                'FoldX_entropy_sidechain', 'FoldX_entropy_mainchain']
    desc += ['Rosetta_dg', 'Rosetta_sasa', 'Rosetta_hbonds']
    # Entropy
    #desc += ['ENM_R6', 'ENM_EXP']
    # Binding
    desc += ['Prodigy_IC_NIS']
    #desc += ['LISA']
    return desc
 

def evaluate(wrkdir, desc_wanted, scores=dict(), force_calc=False):

    desc_set = set()
    for desc in desc_wanted:
        if desc not in scores.keys() or force_calc:
            desc_set.add(desc)
    #print('Wanted : ', desc_wanted)
    #print('Have   : ', scores.keys())
    #print('Calc   : ', desc_set)

    if len(desc_set)>0:
        log.info('Computing new descriptors in %s' % (wrkdir))

    # Molecular descriptors
    if len(desc_set & set(['HB_BH', 'HB_WN', 'HB_KS']))>0:
        scores.update(molecular.hydrogenbond_difference(wrkdir))
    if len(desc_set & set(['BSA', 'BSA_C', 'BSA_A', 'BSA_P', 'NIS_P', 'NIS_C', 'NIS_A', 'NRES']))>0:
        scores.update(molecular.sasa_all(wrkdir))
    if len(desc_set & set(['IC_TOT', 'IC_AA', 'IC_PP', 'IC_CC', 'IC_AP', 'IC_CP', 'IC_AC']))>0:
        scores.update(molecular.intermolecular_contacts(wrkdir))

    # Docking scores
    if 'ZRANK' in desc_set:
        scores.update(docking.zrank(wrkdir))
    if 'ZRANK2' in desc_set:
        scores.update(docking.zrank2(wrkdir))
    if len(desc_set & set(['pyDock', 'pyDock_elec', 'pyDock_vdw', 'pyDock_desolv']))>0:
        scores.update(docking.pydock(wrkdir))
    if 'ATTRACT' in desc_set:
        scores.update(docking.attract(wrkdir))
    if len(desc_set & set(['FireDock', 'FireDock_aVdW', 'FireDock_rVdW', 'FireDock_ACE', 'FireDock_inside',
        'FireDock_aElec', 'FireDock_rElec', 'FireDock_laElec', 'FireDock_lrElec', 'FireDock_hb', 
        'FireDock_piS', 'FireDock_catpiS', 'FireDock_aliph']))>0:
        scores.update(docking.firedock(wrkdir))


    # Statistical Potentials
    if 'RF_HA_SRS' in desc_set:
        scores.update(statistical_potentials.rf_ha_srs(wrkdir))
    if 'RF_CB_SRS_OD' in desc_set:
        scores.update(statistical_potentials.rf_cb_srs_od(wrkdir))
    if 'ipot_aace167' in desc_set:
        scores.update(statistical_potentials.ipot_aace167(wrkdir))
    if 'ipot_aace18' in desc_set:
        scores.update(statistical_potentials.ipot_aace18(wrkdir))
    if 'ipot_aace20' in desc_set:
        scores.update(statistical_potentials.ipot_aace20(wrkdir))
    if 'ipot_rrce20' in desc_set:
        scores.update(statistical_potentials.ipot_rrce20(wrkdir))
    if 'SOAP-PP-Pair' in desc_set:
        scores.update(statistical_potentials.soap_pp(wrkdir))
    if 'SOAP-Protein-OD' in desc_set:
        scores.update(statistical_potentials.soap_protein_od(wrkdir))
    if len(desc_set & set(['DOPE', 'DOPE-HR']))>0:
        scores.update(statistical_potentials.dope(wrkdir))

    # Implicit Solvents
    if 'AGBNP' in desc_set:
        scores.update(implicit_solvents.agbnp(wrkdir))
    if len(desc_set & set(['FACTS_ELEC', 'FACTS_VDW', 'FACTS_GB', 'FACTS_ASP']))>0:
        scores.update(implicit_solvents.facts(wrkdir))
    if len(desc_set & set(['GBMV_ELEC', 'GBMV_VDW', 'GBMV_GB', 'GBMV_ASP']))>0:
        scores.update(implicit_solvents.gbmv(wrkdir))
    if len(desc_set & set(['GBSW_ELEC', 'GBSW_VDW', 'GBSW_GB', 'GBSW_ASP']))>0:
        scores.update(implicit_solvents.gbsw(wrkdir))
    if len(desc_set & set(['CDIE_ELEC', 'CDIE_VDW']))>0:
        scores.update(implicit_solvents.cdie(wrkdir))
    if len(desc_set & set(['RDIE_ELEC', 'RDIE_VDW']))>0:
        scores.update(implicit_solvents.rdie(wrkdir))
    if 'OMM_vacuum' in desc_set:
        scores.update(implicit_solvents.omm_vacuum(wrkdir))
    if 'OMM_HCT' in desc_set:
        scores.update(implicit_solvents.omm_hct(wrkdir))
    if 'OMM_OBC1' in desc_set:
        scores.update(implicit_solvents.omm_obc1(wrkdir))
    if 'OMM_OBC2' in desc_set:
        scores.update(implicit_solvents.omm_obc2(wrkdir))
    if 'OMM_GBn' in desc_set:
        scores.update(implicit_solvents.omm_gbn(wrkdir))
    if 'OMM_GBn2' in desc_set:
        scores.update(implicit_solvents.omm_gbn2(wrkdir))

    # Folding scores
    if len(desc_set & set(['FoldX', 'FoldX_backbone_hbond', 'FoldX_sidechain_hbond', 'FoldX_vdw', 
        'FoldX_elec', 'FoldX_solvation_polar', 'FoldX_solvation_hydrophobic', 'FoldX_entropy_sidechain', 
        'FoldX_entropy_mainchain']))>0:
        scores.update(folding.foldx(wrkdir))
    if len(desc_set & set(['Rosetta_dg', 'Rosetta_sasa', 'Rosetta_hbonds']))>0:
        scores.update(folding.rosetta(wrkdir))

    # Entropy scores
    if 'ENM_R6' in desc_set:
        scores.update(entropy.enm_r6(wrkdir))
    if 'ENM_EXP' in desc_set:
        scores.update(entropy.enm_exp(wrkdir))

    # Binding free energies
    if 'Prodigy_IC_NIS' in desc_set:
        scores.update(binding.prodigy_IC_NIS(wrkdir, desc=scores))
    if 'LISA' in desc_set:
        scores.update(binding.lisa(wrkdir, numthreads=1))

    # Check all descriptors have been calculated
    calculated = set(scores.keys())
    diff = desc_set - calculated
    if len(diff)>0:
        raise ValueError('Unknown descriptors %s. Available descriptors are: %s' % (str(diff), str(all_descriptors())))

    # Return
    return scores

