#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
from openmm.app import *
from openmm import *
from openmm.unit import *
import ppdg
import logging
log = logging.getLogger(__name__)

def get_omm_ener(name, solvent):
    psf    = CharmmPsfFile(name+'-chm.psf')
    pdb    = PDBFile(name+'-chm.pdb')
    params = CharmmParameterSet('/home/simone/opt/ff/charmmff_jul18.prm')
    system = psf.createSystem(
                params,
                nonbondedMethod = CutoffNonPeriodic,
                nonbondedCutoff = 18*angstrom,
                constraints     = None,
                rigidWater      = True,
                hydrogenMass    = 1*amu,
                implicitSolvent = solvent,
                soluteDielectric= 1.0,
                solventDielectric= 80.0
            )
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtosecond)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': "double" , "CudaDeviceIndex" : "0"}
    simulation = Simulation(psf.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=100)
    ener = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie_per_mole)
    return ener

def get_omm_binding(wrkdir, solvent):
    basepath = os.getcwd()
    os.chdir(wrkdir)
    time_start = timer()
    if solvent:
        name = 'OMM_%s' % (solvent)
    else:
        name = 'OMM_vacuum'
    log.info("Getting OpenMM %s scoring..." % (name))
    rec = get_omm_ener('receptor', solvent)
    lig = get_omm_ener('ligand', solvent)
    cpx = get_omm_ener('complex', solvent)
    diff = cpx - lig - rec
    os.chdir(basepath)
    time_end = timer()
    desc = dict()
    desc[name] = diff
    desc['>TIME_%s' % (name)] = time_end - time_start
    return desc

def omm_vacuum(wrkdir):
    return get_omm_binding(wrkdir, None)

def omm_hct(wrkdir):
    return get_omm_binding(wrkdir, HCT)

def omm_obc1(wrkdir):
    return get_omm_binding(wrkdir, OBC1)

def omm_obc2(wrkdir):
    return get_omm_binding(wrkdir, OBC2)

def omm_gbn(wrkdir):
    return get_omm_binding(wrkdir, GBn)

def omm_gbn2(wrkdir):
    return get_omm_binding(wrkdir, GBn2)


if __name__=='__main__':
    ppdg.readconfig('config-ppdg.ini')
    print(get_omm_binding(sys.argv[1], None))
    print(get_omm_binding(sys.argv[1], HCT))
    print(get_omm_binding(sys.argv[1], OBC1))
    print(get_omm_binding(sys.argv[1], OBC2))
    print(get_omm_binding(sys.argv[1], GBn))
    print(get_omm_binding(sys.argv[1], GBn2))

