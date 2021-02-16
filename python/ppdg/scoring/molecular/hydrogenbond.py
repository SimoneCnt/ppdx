#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import mdtraj
import logging
log = logging.getLogger(__name__)

def hydrogenbond(pdbfile):
    """
        Calculate number of hydrogen bonds [1,2] and hydrogen bond energy [3]
        using three methods as implemented in MDtraj [4] for a single PDB file.
        [1] E. N. Baker and R. E. Hubbard, "Hydrogen bonding in globular proteins",
            Progress in Biophysics and Molecular Biology, vol. 44, no. 2, 
            pp. 97-179, 1984.
        [2] P. Wernet et al., "The Structure of the First Coordination Shell in 
            Liquid Water", Science, vol. 304, no. 5673, pp. 995-999, 2004.
        [3] W. Kabsch and C. Sander, "Dictionary of protein secondary structure: 
            Pattern recognition of hydrogen-bonded and geometrical features",
            Biopolymers, vol. 22, no. 12, pp. 2577"2637, 1983.
        [4] R. T. McGibbon et al., "MDTraj: A Modern Open Library for the 
            Analysis of Molecular Dynamics Trajectories", Biophysical Journal, 
            vol. 109, no. 8, pp. 1528-1532, 2015.
    """
    pdb = mdtraj.load(pdbfile)
    desc = dict()
    desc['HB_BH'] = float(mdtraj.baker_hubbard(pdb).shape[0])
    desc['HB_WN'] = float(mdtraj.wernet_nilsson(pdb)[0].shape[0])
    desc['HB_KS'] = float(mdtraj.kabsch_sander(pdb)[0].sum())
    return desc

def hydrogenbond_difference(wrkdir, cpxname='complex.pdb', recname='receptor.pdb', ligname='ligand.pdb'):
    """
        Calculate difference in number of hydrogen bonds [1,2] and hydrogen 
        bond energy [3] using three methods as implemented in MDtraj [4] for a 
        protein-protein complex upon binding.
        [1] E. N. Baker and R. E. Hubbard, "Hydrogen bonding in globular proteins",
            Progress in Biophysics and Molecular Biology, vol. 44, no. 2, 
            pp. 97-179, 1984.
        [2] P. Wernet et al., "The Structure of the First Coordination Shell in 
            Liquid Water", Science, vol. 304, no. 5673, pp. 995-999, 2004.
        [3] W. Kabsch and C. Sander, "Dictionary of protein secondary structure: 
            Pattern recognition of hydrogen-bonded and geometrical features",
            Biopolymers, vol. 22, no. 12, pp. 2577"2637, 1983.
        [4] R. T. McGibbon et al., "MDTraj: A Modern Open Library for the 
            Analysis of Molecular Dynamics Trajectories", Biophysical Journal, 
            vol. 109, no. 8, pp. 1528-1532, 2015.
    """
    time_start = timer()
    log.info("Getting hydrogen bonds: HB_BH, HB_WN and HB_KS")
    cpx = hydrogenbond(os.path.join(wrkdir, cpxname))
    lig = hydrogenbond(os.path.join(wrkdir, ligname))
    rec = hydrogenbond(os.path.join(wrkdir, recname))
    desc = dict()
    for key in cpx.keys():
        desc[key] = cpx[key] - rec[key] - lig[key]
    time_end = timer()
    desc['>TIME_hydrogenbond'] = time_end-time_start
    return desc


