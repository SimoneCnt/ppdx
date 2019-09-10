#!/usr/bin/env python3

import sys, os
from timeit import default_timer as timer
import numpy as np
import mdtraj
import logging
log = logging.getLogger(__name__)

LEVY_CLASS_INTERIOR = 1
LEVY_CLASS_SURFACE = 2
LEVY_CLASS_SUPPORT = 3
LEVY_CLASS_RIM = 4
LEVY_CLASS_CORE = 5

def sasa_mdtraj(pdbfile):
    """
        Calculate the SASA per residue using the Shrake & Rupley algorithm [1]
        as implemented in mdtraj [2]. Return values in angstrom^2.
        [1] A. Shrake and J. A. Rupley, "Environment and exposure to solvent 
            of protein atoms. Lysozyme and insulin", Journal of Molecular 
            Biology, vol. 79, no. 2, pp. 351-371, 1973.
        [2] R. T. McGibbon et al., "MDTraj: A Modern Open Library for the 
            Analysis of Molecular Dynamics Trajectories", Biophysical Journal, 
            vol. 109, no. 8, pp. 1528-1532, 2015.
    """
    pdb = mdtraj.load(pdbfile)
    res = [ residue.name for residue in pdb.topology.residues ]
    sasa = mdtraj.shrake_rupley(pdb, mode='residue')[0]*100.0
    return res, sasa

def sasa(pdbfile):
    """
        Calculate the SASA per residue. Return values in angstrom^2.
        For the moment, just an alias to sasa_mdtraj.
    """
    return sasa_mdtraj(pdbfile)

def sasa_relative(pdbfile, method='miller'):
    """
        Calculate the relative SASA per residue. The reference SASA are taken 
        from Gly-X-Gly tripeptide as calculated in S.Miller [1], or 
        C.Chothia [2].
        [1] S. Miller, J. Janin, A. M. Lesk, and C. Chothia, "Interior and 
            surface of monomeric proteins", Journal of Molecular Biology, 
            vol. 196, no. 3, pp. 641-656, 1987.
        [2] C. Chothia, "The nature of the accessible and buried surfaces in 
            proteins", Journal of Molecular Biology, vol. 105, no. 1, 
            pp. 1-12, 1976.
    """
    if method=='miller':
        ref = {
            'ALA':113, 'ARG':241, 'ASN':158, 'ASP':151, 'CYS':140,
            'GLN':189, 'GLU':183, 'GLY': 85, 'HIS':194, 'ILE':182,
            'LEU':180, 'LYS':211, 'MET':204, 'PHE':218, 'PRO':143,
            'SER':122, 'THR':146, 'TRP':259, 'TYR':229, 'VAL':160
            }
    elif method=='chothia':
        ref = {
            'ALA':115, 'ARG':225, 'ASN':160, 'ASP':150, 'CYS':135,
            'GLN':180, 'GLU':190, 'GLY': 75, 'HIS':195, 'ILE':175,
            'LEU':170, 'LYS':200, 'MET':185, 'PHE':210, 'PRO':145,
            'SER':115, 'THR':140, 'TRP':255, 'TYR':230, 'VAL':155
            }
    else:
        raise ValueError('Unknown method <%s> to calculate the relative SASA. Use miller or chothia' % (method))
    res, asa = sasa(pdbfile)
    rel = [ s/ref[r] for r,s in zip(res, asa) ]
    return res, asa, rel

def sasa_difference(wrkdir, cpxname='complex.pdb', recname='receptor.pdb', ligname='ligand.pdb'):
    """
        Calculate the difference in SASA upon binding.
    """
    res_cpx, sasa_cpx = sasa(os.path.join(wrkdir, cpxname))
    res_lig, sasa_lig = sasa(os.path.join(wrkdir, ligname))
    res_rec, sasa_rec = sasa(os.path.join(wrkdir, recname))
    sasa_bound = np.array(sasa_cpx)
    sasa_unbound = np.concatenate((sasa_rec, sasa_lig), axis=None)
    sasa_diff = sasa_unbound - sasa_bound
    return res_cpx, sasa_bound, sasa_unbound, sasa_diff

def sasa_relative_difference(wrkdir, cpxname='complex.pdb', recname='receptor.pdb', ligname='ligand.pdb'):
    """
        Calculate the difference in relative SASA upon binding.
    """
    res_cpx, sasa_cpx, rsasa_cpx = sasa_relative(os.path.join(wrkdir, cpxname))
    res_lig, sasa_lig, rsasa_lig = sasa_relative(os.path.join(wrkdir, ligname))
    res_rec, sasa_rec, rsasa_rec = sasa_relative(os.path.join(wrkdir, recname))
    sasa_bound = np.array(sasa_cpx)
    sasa_unbound = np.concatenate((sasa_rec, sasa_lig), axis=None)
    sasa_diff = sasa_unbound - sasa_bound
    rsasa_bound = np.array(rsasa_cpx)
    rsasa_unbound = np.concatenate((rsasa_rec, rsasa_lig), axis=None)
    rsasa_diff = rsasa_unbound - rsasa_bound
    return res_cpx, sasa_bound, sasa_unbound, sasa_diff, rsasa_bound, rsasa_unbound, rsasa_diff
 
def sasa_all(wrkdir, cpxname='complex.pdb', recname='receptor.pdb', ligname='ligand.pdb'):
    """
        Calculate the difference in SASA upon binding aka Buried Surface Area (BSA)
        as total, polar, apolar and charged, and the Percentage of Non-Interacting 
        Surface (NIS) as polar, apolar and charged.
        Return everything as dictionary
        For more info on NIS see [1].
        [1] P. L. Kastritis, J. P. G. L. M. Rodrigues, G. E. Folkers, R. Boelens, 
            and A. M. J. J. Bonvin, "Proteins Feel More Than They See: 
            Fine-Tuning of Binding Affinity by Properties of the Non-Interacting 
            Surface," Journal of Molecular Biology, vol. 426, no. 14, pp. 2632–2652, 2014.
    """
    time_start = timer()
    log.info("Getting SASA, BSA and NIS")
    aaprop = {
        'ALA':'A', 'CYS':'P', 'GLU':'C', 'ASP':'C', 'GLY':'A',
        'PHE':'A', 'ILE':'A', 'HIS':'P', 'LYS':'C', 'MET':'A',
        'LEU':'A', 'ASN':'P', 'GLN':'P', 'PRO':'A', 'SER':'P',
        'ARG':'C', 'THR':'P', 'TRP':'P', 'VAL':'A', 'TYR':'P'
        }
    #atoms_per_residue = {
    #    'ALA':1, 'CYS':2, 'GLU':5, 'ASP':4, 'GLY':0,
    #    'PHE':7, 'ILE':4, 'HIS':6, 'LYS':5, 'MET':4,
    #    'LEU':4, 'ASN':4, 'GLN':5, 'PRO':4, 'SER':2,
    #    'ARG':7, 'THR':3, 'TRP':10, 'VAL':3, 'TYR':8
    #    }
    res_cpx, sasa_bound, sasa_unbound, sasa_diff, rsasa_bound, rsasa_unbound, rsasa_diff = sasa_relative_difference(wrkdir, cpxname, recname, ligname)
    desc = {
        'BSA'   : 0.0, # 'BURIED_ATOMS': 0.0,
        'BSA_P' : 0.0, 'BSA_A' : 0.0, 'BSA_C' : 0.0,
        'NIS_P' : 0.0, 'NIS_A' : 0.0, 'NIS_C' : 0.0
        }
    for res, asa in zip(res_cpx, sasa_diff):
        k = 'BSA_'+aaprop[res]
        desc[k] += asa
        desc['BSA'] += asa
        #desc['BURIED_ATOMS'] += atoms_per_residue[res]
    for res, asa in zip(res_cpx, rsasa_bound):
        if asa>0.05:
            k = 'NIS_'+aaprop[res]
            desc[k] += 1
    total = desc['NIS_P'] + desc['NIS_A'] + desc['NIS_C']
    desc['NIS_P'] *= 100.0/total
    desc['NIS_A'] *= 100.0/total
    desc['NIS_C'] *= 100.0/total
    desc['NRES'] = len(res_cpx) # Number of residues in the complex
    time_end = timer()
    desc['>TIME_sasa'] = time_end - time_start
    return desc
    
def classify_residues_levy(wrkdir, cpxname='complex.pdb', recname='receptor.pdb', ligname='ligand.pdb'):
    """
        Classify the residues of the protein-protein complex in interior,
        surface, support, rim, or core according to Levy definition [1].
        [1] E. D. Levy, "A Simple Definition of Structural Regions in Proteins 
            and Its Use in Analyzing Interface Evolution", Journal of 
            Molecular Biology, vol. 403, no. 4, pp. 660-670, 2010.
    """
    res_cpx, _, _, _, rsasa_bound, rsasa_unbound, rsasa_diff = sasa_relative_difference(wrkdir, cpxname, recname, ligname)
    klass = list()
    for b, u, d in zip(rsasa_bound, rsasa_unbound, rsasa_diff):
        if d<0.01:
            if b<0.25:
                klass.append(LEVY_CLASS_INTERIOR)
            else:
                klass.append(LEVY_CLASS_SURFACE)
        else:
            if u>0.25:
                if b>0.25:
                    klass.append(LEVY_CLASS_RIM)
                else:
                    klass.append(LEVY_CLASS_CORE)
            else:
                if b>0.25:
                    log.error('This condition is impossible! Bound rSASA=%.3f cannot be greater than unbound rSASA=%.3f!' % (b, u))
                    klass.append(None)
                else:
                    klass.append(LEVY_CLASS_SUPPORT)
    return res_cpx, klass


def stickiness(wrkdir):
    """
        Compute the average stickiness of the core interface residues according
        to Levy [1].
        [1] E. D. Levy, S. De, and S. A. Teichmann, "Cellular crowding imposes 
            global constraints on the chemistry and evolution of proteomes", 
            PNAS, vol. 109, no. 50, pp. 20461-20466, 2012.
    """
    stk = {
        'ALA': 0.0062, 'CYS': 1.0372, 'ASP':-0.7485, 'GLU':-0.7893, 'PHE': 1.2727,
        'GLY':-0.1771, 'HIS': 0.1204, 'ILE': 1.1109, 'LYS':-1.1806, 'LEU': 0.9138,
        'MET': 1.0124, 'ASN':-0.2693, 'PRO':-0.1799, 'GLN':-0.4114, 'ARG':-0.0876,
        'SER': 0.1376, 'THR': 0.1031, 'VAL': 0.7599, 'TRP': 0.7925, 'TYR': 0.8806
        }
    time_start = timer()
    log.info("Getting Stickiness")
    res, klass = classify_residues_levy(wrkdir)
    stks = [ stk[r] for r, k in zip(res, klass) if k==LEVY_CLASS_CORE ]
    desc = dict()
    desc['sticky_tot'] = np.sum(stks)
    desc['sticky_avg'] = np.average(stks)
    time_end = timer()
    desc['>TIME_stickiness'] = time_end - time_start
    return desc


if __name__=='__main__':
    if len(sys.argv)==2:
        print(sasa_all(sys.argv[1]))
        print(stickiness(sys.argv[1]))
    elif len(sys.argv)==5:
        print(sasa_all(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
        print(stickiness(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))
    else:
        print('Usage: %s wrkdir [complex.pdb receptor.pdb ligand.pdb]' % (sys.argv[0]))
        quit()


# To check AAindex
# https://www.genome.jp/aaindex/AAindex/list_of_indices


