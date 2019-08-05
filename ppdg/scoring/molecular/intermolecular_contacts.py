#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import Bio.PDB
import logging
log = logging.getLogger(__name__)

def intermolecular_contacts(wrkdir, cpxname='complexAB.pdb', d_cutoff=5.5):
    """
        Calculates intermolecular contacts between the two chains in a PDB file
        based on the definition in Vangone & Bonvin [1]. The residues are then
        classified by residue polarity (polar, apolar, charged; six class in 
        total). Inspired by the predict_ic.py code of Vangone & Bonvin.
        WARN! The code has been adapted to count the IC contacts just between
            two chains of a protein-protein complex, thus the PDB file has to 
            contain exactly 2 chains. If ti contains more, but still one set is
            the receptor and the other the ligand, just renaim the chains so
            that the full receptor is chain A and the full ligand is chain B.
        [1] A. Vangone and A. M. Bonvin, "Contacts-based prediction of binding 
            affinity in protein–protein complexes", eLife, vol. 4, p. e07454, 2015.
    """
    time_start = timer()
    log.info('Getting intermolecular contacts (IC)')
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure('complex', os.path.join(wrkdir, cpxname))
    atom_list = list(structure.get_atoms())
    _ignore = lambda x: x.element == 'H'
    for atom in atom_list:
        if _ignore(atom):
            residue = atom.parent
            residue.detach_child(atom.name)
    atom_list = list(structure.get_atoms())
    ns = Bio.PDB.NeighborSearch(atom_list)
    all_list = ns.search_all(radius=d_cutoff, level='R')
    ic_list = [c for c in all_list if c[0].parent.id != c[1].parent.id]
    aa_character_ic = {
        'ALA':'A', 'CYS':'A', 'GLU':'C', 'ASP':'C', 'GLY':'A',
        'PHE':'A', 'ILE':'A', 'HIS':'C', 'LYS':'C', 'MET':'A',
        'LEU':'A', 'ASN':'P', 'GLN':'P', 'PRO':'A', 'SER':'P',
        'ARG':'C', 'THR':'P', 'TRP':'A', 'VAL':'A', 'TYR':'A'
        }
    bins = { 
        'IC_TOT': len(ic_list),
        'IC_AA': 0, 'IC_PP': 0,
        'IC_CC': 0, 'IC_AP': 0,
        'IC_CP': 0, 'IC_AC': 0,
        }
    for (res_i, res_j) in ic_list:
        contact_type = 'IC_'+''.join(sorted((aa_character_ic.get(res_i.resname), aa_character_ic.get(res_j.resname))))
        bins[contact_type] += 1
    time_end = timer()
    bins['>TIME_intermolecular_contacts'] = time_end - time_start
    return bins

