#!/usr/bin/env python3

import os, sys
import numpy as np
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def pydock(wrkdir):
    """
        Calculate pyDock scoring.
        [1] T. M.-K. Cheng, T. L. Blundell, and J. Fernandez-Recio, 
            "pyDock: Electrostatics and desolvation for effective scoring of 
            rigid-body protein–protein docking", Proteins: Structure, 
            Function, and Bioinformatics, vol. 68, no. 2, pp. 503-515, 2007.
    """
    time_start = timer()
    log.info("Getting pyDock scoring...")
    basepath = os.getcwd()
    os.chdir(wrkdir)
    rec = ppdg.Pdb('receptor.pdb')
    rec.set_chain('A')
    rec.write('receptorA.pdb')
    lig = ppdg.Pdb('ligand.pdb')
    lig.set_chain('B')
    lig.write('ligandB.pdb')
    xavg = np.mean([ atom.x for atom in lig.atoms])
    yavg = np.mean([ atom.y for atom in lig.atoms])
    zavg = np.mean([ atom.z for atom in lig.atoms])
    with open("pydock.ini", 'w') as fp:
        fp.write("[receptor]\n")
        fp.write("pdb     = receptorA.pdb\n")
        fp.write("mol     = %s\n" % ('A'))
        fp.write("newmol  = %s\n" % ('A'))
        fp.write("\n")
        fp.write("[ligand]\n")
        fp.write("pdb     = ligandB.pdb\n")
        fp.write("mol     = %s\n" % ('B'))
        fp.write("newmol  = %s\n" % ('B'))
    stdout, stderr, ret = ppdg.tools.execute("%s pydock setup" % (os.path.join(ppdg.PYDOCK, 'pyDock3')))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("pyDock setup failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))
    with open("pydock.rot", 'w') as fp:
        fp.write("1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 %f %f %f 1\n" % (xavg, yavg, zavg))
    stdout, stderr, ret = ppdg.tools.execute("%s pydock dockser" % (os.path.join(ppdg.PYDOCK, 'pyDock3')))
    if ret!=0:
        os.chdir(basepath)
        raise ValueError("pyDock dockser failed! Returned code is %d\nSTDOUT:\n%s\nSTDERR:\n%s" % (ret, stdout, stderr))
    with open('pydock.ene', 'r') as fp:
        line = fp.readlines()[-1]
        _, elec, desolv, vdw, total, _ = line.split()
    os.chdir(basepath)
    desc = {'pyDock':float(total), 'pyDock_elec':float(elec), 'pyDock_desolv':float(desolv), 'pyDock_vdw':float(vdw)}
    time_end = timer()
    desc['>TIME_pyDock'] = time_end - time_start
    return desc

if __name__=='__main__':
    ppdg.readconfig('config-ppdg.ini')
    print(pydock(sys.argv[1]))

