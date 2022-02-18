#!/usr/bin/env python3

from timeit import default_timer as timer
import Bio.PDB
import logging
log = logging.getLogger(__name__)
from ppdx.scoring.molecular import sasa_all, intermolecular_contacts

def prodigy_IC_NIS(wrkdir, cpxname='complexAB.pdb', recname='receptor.pdb', ligname='ligand.pdb', desc=dict()):
    """
        Calculates the predicted binding affinity value
        based on the IC-NIS model.
        https://github.com/haddocking/binding_affinity
    """
    log.info('Getting Prodigy_IC_NIS')
    desc_set = set(desc.keys())
    if len(set(['IC_CC', 'IC_AC', 'IC_PP', 'IC_AP'])-desc_set)>0:
        desc.update(intermolecular_contacts(wrkdir, cpxname))
    if len(set(['NIS_A', 'NIS_C'])-desc_set)>0:
        desc.update(sasa_all(wrkdir, cpxname, recname, ligname))
    time_start = timer()
    ic_cc = desc['IC_CC']
    ic_ca = desc['IC_AC']
    ic_pp = desc['IC_PP']
    ic_pa = desc['IC_AP']
    nis_a = desc['NIS_A']
    nis_c = desc['NIS_C']
    desc['Prodigy_IC_NIS'] = -0.09459*ic_cc + -0.10007*ic_ca + 0.19577*ic_pp + -0.22671*ic_pa \
                                + 0.18681*nis_a + 0.13810*nis_c + -15.9433
    time_end = timer()
    desc['>TIME_prodigy_ic_nis'] = time_end - time_start
    return desc

#def prodigy_NIS(wrkdir, cpxname='complexAB.pdb', recname='receptor.pdb', ligname='ligand.pdb', desc=dict()):
#    """
#        Calculates the predicted binding affinity value
#        based on the NIS model.
#    """
#    desc_set = set(desc.keys())
#    if len(set(['NIS_P', 'NIS_C', 'BURIED_ATOMS'])-desc_set)>0:
#        desc.update(sasa_all(wrkdir, cpxname, recname, ligname))
#    nis_c = desc['NIS_C']
#    nis_p = desc['NIS_P']
#    natoms = desc['BURIED_ATOMS']
#    desc['Prodigy_NIS'] = 0.0856851248873*nis_p + -0.0685254498746*nis_c + 0.0261591389985*natoms + 3.0124939659498
#    return desc

if __name__=='__main__':
    import sys, os
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if len(sys.argv)==2:
        desc = prodigy_IC_NIS(sys.argv[1])
        desc = prodigy_NIS(sys.argv[1], desc=desc)
    elif len(sys.argv)==5:
        desc = prodigy_IC_NIS(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        desc = prodigy_NIS(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], desc=desc)
    else:
        print('Usage: %s wrkdir [complex.pdb receptor.pdb ligand.pdb]' % (sys.argv[0]))
        quit()
    print(desc)

