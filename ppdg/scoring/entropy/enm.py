#!/usr/bin/env python3

import os
from timeit import default_timer as timer
import ppdg
from math import sqrt, exp
import numpy as np
import logging
log = logging.getLogger(__name__)

def enm_exp(wrkdir):
    dS, time = enm_dS(wrkdir, 'exponential')
    desc = dict()
    desc['ENM_EXP'] = dS
    desc['>TIME_ENM_EXP'] = time
    return desc

def enm_r6(wrkdir):
    dS, time = enm_dS(wrkdir, 'overR6')
    desc = dict()
    desc['ENM_R6'] = dS
    desc['>TIME_ENM_R6'] = time
    return desc

def enm_dS(wrkdir, method):
    time_start = timer()
    if not os.path.isfile(os.path.join(wrkdir, 'complex.pdb')):
        raise ValueError('File complex.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'ligand.pdb')):
        raise ValueError('File ligand.pdb does not exist in %s.' % (wrkdir))
    if not os.path.isfile(os.path.join(wrkdir, 'receptor.pdb')):
        raise ValueError('File receptor.pdb does not exist in %s.' % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting ENM %s entropy..." % (method))
    cpx = enm_entropy('complex.pdb', K=1000.0, cutoff=999.0, spring=method)
    lig = enm_entropy('ligand.pdb', K=1000.0, cutoff=999.0, spring=method)
    rec = enm_entropy('receptor.pdb', K=1000.0, cutoff=999.0, spring=method)
    os.chdir(basepath)
    time_end = timer()
    return cpx-lig-rec, time_end - time_start

def enm_entropy(fname, K=50.0, cutoff=7.0, spring='constant'):
    """
        Given a pdb file, make the Elastic Network Model (ENM), 
        build its hessian matrix, diagonalize it, and return its entropy.
    """

    pdb_all = ppdg.Pdb(fname)
    pdb_all.remove_hydrogens()
    pdb = [ [atom.x, atom.y, atom.z] for atom in pdb_all.atoms ]

    def spring_constant(k,d):
        return k/(d*d)
    def spring_exponential(k,d,d0=7.0):
        return k*exp(-(d/d0)**2)
    def spring_overR6(k,d):
        return k/(d**4)

    if spring=='constant':
        spr = spring_constant
    elif spring=='exponential':
        spr = spring_exponential
    elif spring=='overR6':
        spr = spring_overR6
    else:
        log.error("Unkown spring method <%s>. Available are: constant, exponential, overR6" % spring)
        return float('nan')

    #print("Making ENM")
    enm = []
    natoms = len(pdb)
    for i in range(natoms):
        for j in range(i+1, natoms):
            dx = pdb[i][0] - pdb[j][0]
            dy = pdb[i][1] - pdb[j][1]
            dz = pdb[i][2] - pdb[j][2]
            d  = sqrt(dx*dx+dy*dy+dz*dz)
            if (d<cutoff):
                enm.append([i, j, dx, dy, dz, d])
    #print("Among all N*(N-1)/2 = %d pairs, %d were selected." % (natoms*(natoms-1)/2, len(enm)))
    #print("Building Hessian Matrix")
    nfree = natoms*3
    hess = np.zeros((nfree, nfree))
    for bond in enm:
        i  = bond[0]
        j  = bond[1]
        dx = bond[2]
        dy = bond[3]
        dz = bond[4]
        d  = bond[5]
        tmp = spr(K,d)
        for m,v1 in enumerate([dx,dy,dz]):
            for n,v2 in enumerate([dx,dy,dz]):
                hess[3*i+m][3*i+n] += tmp*v1*v2
                hess[3*j+m][3*j+n] += tmp*v1*v2
                hess[3*i+m][3*j+n] -= tmp*v1*v2
                hess[3*j+m][3*i+n] -= tmp*v1*v2
    #print("Diagonalizing it...")
    #eival, eivec = np.linalg.eig(hess)
    eival = np.linalg.eigvalsh(hess)
    #print("Calculating the entropy")
    S = 0
    skipped = 0
    toout=""
    for i in eival:
        if i>1E-6:
            S+= 1 + np.log(sqrt(i))
        else:
            skipped+=1
            toout += "%d %f\n" % (skipped, i)
    if skipped!=6:
        log.warning("Warning! There were %d skipped frequencees below 1E-4! Exactly 6 were expected!" % skipped)
        log.warning(toout)
    return S

