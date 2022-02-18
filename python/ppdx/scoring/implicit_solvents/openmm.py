#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import logging
log = logging.getLogger(__name__)

def get_omm_ener(name, solvent):
    import openmm
    import openmm.app
    import openmm.unit as unit
    if solvent=='HCT': solvent=openmm.app.HCT
    if solvent=='OBC1': solvent=openmm.app.OBC1
    if solvent=='OBC2': solvent=openmm.app.OBC2
    if solvent=='GBn': solvent=openmm.app.GBn
    if solvent=='GBn2': solvent=openmm.app.GBn2
    psf    = openmm.app.CharmmPsfFile(name+'-chm.psf')
    pdb    = openmm.app.PDBFile(name+'-chm.pdb')
    params = openmm.app.CharmmParameterSet('/home/simone/opt/ff/charmmff_jul18.prm')
    system = psf.createSystem(
                params,
                nonbondedMethod = openmm.app.CutoffNonPeriodic,
                nonbondedCutoff = 18*unit.angstrom,
                constraints     = None,
                rigidWater      = True,
                hydrogenMass    = 1*unit.amu,
                implicitSolvent = solvent,
                soluteDielectric= 1.0,
                solventDielectric= 80.0
            )
    integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 2*unit.femtosecond)
    platform = openmm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': "double" , "CudaDeviceIndex" : "0"}
    simulation = openmm.app.Simulation(psf.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(maxIterations=100)
    ener = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
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
    return get_omm_binding(wrkdir, 'HCT')

def omm_obc1(wrkdir):
    return get_omm_binding(wrkdir, 'OBC1')

def omm_obc2(wrkdir):
    return get_omm_binding(wrkdir, 'OBC2')

def omm_gbn(wrkdir):
    return get_omm_binding(wrkdir, 'GBn')

def omm_gbn2(wrkdir):
    return get_omm_binding(wrkdir, 'GBn2')


