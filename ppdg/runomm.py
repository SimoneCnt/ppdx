#!/usr/bin/env python3

import sys
import os
import simtk.openmm.app as app
import simtk.openmm as omm
import simtk.unit as u
import openmmscripts as omms
from math import sqrt
import numpy as np

import logging
log = logging.getLogger(__name__)

def runomm(
        # Starting topology
        psffile             = None,
        # Starting coordinates (either pdb or cor)
        pdbfile             = None,
        crdfile             = None,
        # Force field
        prmfile             = None,
        toppar              = None,
        # Simulation timestep and total simulation time
        timestep            = 2*u.femtoseconds,
        simultime           = 0*u.picoseconds,
        # Time interval for writing outputs (dcd, log, ene)
        outtime             = 10*u.picoseconds,
        # Base name for output files
        outname             = None,
        # Filename of a xml restart file
        prevxml             = None,
        # System properties
        temperature         = 300*u.kelvin,
        pressure            = 1*u.atmosphere,
        pressure_frequency  = 25,
        cutoffmethod        = app.PME,
        cutoff              = 14*u.angstrom,
        switchdist          = None,
        constraints         = app.HBonds,
        hydrogenmass        = None,
        boxx                = None,
        boxy                = None,
        boxz                = None,
        implicit_solv       = None,
        dielec_solute       = 1.0,
        dielec_solvent      = 80.0,
        implicit_kappa      = None,
        platform            = 'CUDA',
        CUDAprec            = 'mixed',
        CUDAid              = '0',
        # Energy minimization. If prevxml is define, minimization is skipped
        domini              = True,
        ministep            = 1000,
        minitol             = 10*u.kilojoule/u.mole,
        # Do single point energy
        singlepoint         = False,
        absolute_force      = None,
        plumedscript        = None,
        smt                 = False,
        adt                 = False,
        variable_dt         = False,
        # Alchemical Thermodynamic Integration and Dynamics
        alch                = False,
        alch_pdb            = None,
        alch_lambda         = 1.0,
        alch_dlambda        = 0.005,
        alch_outfile        = 'alchfile.ti',
        alch_outtime        = 1*u.picoseconds,
        alch_dynamics       = False,
        alch_dt             = 1,
        alch_beta           = 1,
        # Dynamo plugin: String method (smcv) and water shell (wshell)
        dynamo_inp          = 'dynamo.inp',
        dynamo_out          = 'dynamo.out',
        dynamo_pdb          = None,
        smcv                = None,
        wshell              = None,
        # Fixed Sphere
        fixsphere           = False,
        fixsphere_center    = None,
        fixsphere_radius    = None,
        # Implicit Surface
        isurf               = False,
        isurf_de            = 4*u.kilojoule_per_mole,
        isurf_a             = 2*(1/u.angstrom),
        isurf_z0            = 5*u.angstrom,
        isurf_zmax          = 100*u.angstrom,
        # REST2
        rest2_thot          = None,
        rest2_atoms         = None,
        **kwargs
    ):

    # Check for unrecognized options
    for key, val in kwargs.items():
        log.warning("Unknown option %s = %s" % (key, val))

    # Set the output file names
    if not outname:
        raise ValueError("You must specify and outname!")
    ommlogfile = outname+'.log'
    dcdoutfile = outname+'.dcd'
    eneoutfile = outname+'.ene'
    pdboutfile = outname+'-final.pdb'
    xmloutfile = outname+'-final.xml'

    # Read topology
    if not psffile:
        log.error("Did not set a topology file! This is needed! Please pass it using the <psffile> variable.")
        return None
    log.info("Reading topology file %s" % psffile)
    psf = app.CharmmPsfFile(psffile)

    # Read coordinated from charmm crd or from pdb
    if pdbfile:
        log.info("Reading PDB %s " % pdbfile)
        pdb = app.PDBFile(pdbfile)
    elif crdfile:
        log.info("Reading CHARMM CRD %s " % crdfile)
        pdb = app.CharmmCrdFile(crdfile)
    else:
        log.error("Neither pdbfile nor crdfile is set! Use <pdbfile> or <crdfile> variables.")
        return None

    # Get box size
    if not (boxx and boxy and boxz):
        log.info("No box size found. Guessing it from coordinates")
        coor = np.array(pdb.positions.value_in_unit(u.nanometer))
        boxx = boxy = boxz = np.amax(coor) - np.amin(coor) + 0.1
    log.info("Setting box size: %s %s %s" % (boxx, boxy, boxz))
    psf.setBox(boxx, boxy, boxz)

    # Read Force field parameters read from $HOME/opt/ff/
    if prmfile:
        log.info("Readinf parameter file %s" % (prmfile))
        params = app.CharmmParameterSet(prmfile)
    else:
        log.info("Reading CHARMM36 parameter file")
        if toppar:
            params = app.CharmmParameterSet(*toppar)
        else:
            params = app.CharmmParameterSet(os.environ['HOME']+'/opt/ff/merged.prm')

    # Simulation parameters
    log.info("Creating the system!")
    log.info("Cutoff: cutoffmethod " + str(cutoffmethod) + " cutoff " + str(cutoff) + " switchdist " + str(switchdist))
    log.info("Implicit solvent " + str(implicit_solv) + " implicitKappa " + str(implicit_kappa) + " dielectric " + str(dielec_solute) + " " + str(dielec_solvent))
    log.info("Constraints " + str(constraints) + " Hydrogenmass " + str(hydrogenmass))
    system = psf.createSystem(
                    params, 
                    nonbondedMethod = cutoffmethod,
                    nonbondedCutoff = cutoff, 
                    switchDistance  = switchdist,
                    constraints     = constraints,
                    rigidWater      = True,
                    hydrogenMass    = hydrogenmass,
                    removeCMMotion  = False,
                    implicitSolvent = implicit_solv, 
                    soluteDielectric= dielec_solute,
                    solventDielectric= dielec_solvent,
                    implicitSolventKappa = implicit_kappa,
                )

    # Integrator
    if (variable_dt==False):
        log.info("Setup Langevin Integrator with timestep " + str(timestep))
        integrator = omm.LangevinIntegrator(temperature, 10/u.picosecond, timestep)
    else:
        log.info("Setup Langevin Integrator with variable timestep")
        integrator = omm.VariableLangevinIntegrator(temperature, 10/u.picosecond, 0.001)

    # Barostat
    if pressure:
        log.info("Adding Monte Carlo Barostat with pressure %s and frequencey %d" % (str(pressure), pressure_frequency))
        barostat = omm.MonteCarloBarostat(pressure, temperature, pressure_frequency)
        system.addForce(barostat)

    # REST2
    if rest2_thot and rest2_atoms:
        log.info("Using REST2 with low temperature of %s and hot temperature of %s on a set of %d atoms: %s" % (str(temperature), str(rest2_thot), len(rest2_atoms), str(rest2_atoms)))
        omms.rest2(system, temperature, rest2_thot, rest2_atoms)

    # Absolute harmonic force on CA ~ on corfile positions
    if absolute_force:
        log.info("Adding absolute harmonic force on CA with k = " + str(absolute_force))
        system.addForce(omms.AbsoluteForce(pdb, absolute_force))

    # Fix Sphere: fix in position and set non-bonded parameters to zero for all atoms at distance higher than fixsphere_radius from fixsphere_center (atom id)
    if fixsphere:
        omms.fixsphere(system, pdb, fixsphere_center, fixsphere_radius)

    # Plumed 
    if plumedscript:
        log.info("Adding plumed force")
        from openmmplumed import PlumedForce
        plumedforce = PlumedForce(plumedscript)
        plumedforce.setForceGroup(omms.BIAS_FORCE_GROUP)
        system.addForce(plumedforce)

    # Implicit surface
    if isurf:
        isurf = omm.CustomExternalForce("De*(exp(-2*a*(z-z0))-2*exp(-a*(z-z0)))-10/(z-zmax)")
        isurf.addGlobalParameter("De", isurf_de.value_in_unit(u.kilojoule_per_mole))
        isurf.addGlobalParameter("a", isurf_a.value_in_unit((u.nanometer)**(-1)))
        isurf.addGlobalParameter("z0", isurf_z0.value_in_unit(u.nanometer))
        isurf.addGlobalParameter("zmax", isurf_zmax.value_in_unit(u.nanometer))
        isurf.setForceGroup(omms.BIAS_FORCE_GROUP)
        for i in range(system.getNumParticles()):
            isurf.addParticle(i)
        system.addForce(isurf)


    # Dynamo plugin for water shell and string
    if smcv or wshell:
        log.info("Adding Dynamo Force")
        from openmmdynamo import DynamoForce
        omms.dynamo_writeinp(dynamo_inp, dynamo_pdb, smcv, wshell)
        system.addForce(DynamoForce(dynamo_inp, dynamo_out))

    # Platform
    pltform = omm.Platform.getPlatformByName(platform)
    properties = {}
    if platform=='CUDA':
        properties = {'CudaPrecision': CUDAprec , "CudaDeviceIndex" : str(CUDAid)}

    # Simulation
    if alch:
        log.info("Preparing simulation for alchemical system")
        alch_system = omms.make_alchemical_system(system, alch_pdb, alch_lambda, annihilate=False)
        simulation = app.Simulation(psf.topology, alch_system, integrator, pltform, properties)
    else:
        log.info("Preparing simulation")
        simulation = app.Simulation(psf.topology, system, integrator, pltform, properties)

    # Set coordinates and velocities
    if prevxml:
        log.info("Reading coordinates and velocities from %s " % prevxml)
        #with open(prevxml, 'r') as f:
        #    xml=f.read();
        #    oldstate=omm.XmlSerializer.deserialize(xml)
        #    simulation.context.setPositions(oldstate.getPositions());
        #    simulation.context.setVelocities(oldstate.getVelocities());
        simulation.loadState(prevxml)
    else:
        log.info("Using initial coordinates and random velocities.")
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(temperature)

    # Do energy minimization
    if domini and not prevxml:
        energy = omms.get_energy(simulation.context)
        log.info('Energy Before Minimization is %s' % (str(energy)))
        omms.print_energy_decomposition(simulation.context)
        log.info("Running energy minimization (maxsteps = %s tolerance = %s)" % (ministep, minitol))
        simulation.minimizeEnergy(maxIterations=ministep, tolerance=minitol)
        energy = omms.get_energy(simulation.context)
        log.info('Energy After Minimization is %s' % energy)
        omms.print_energy_decomposition(simulation.context)

    # Single point energy
    if singlepoint:
        energy = omms.get_energy(simulation.context)
        log.info('Single point energy is %s' % (str(energy)))
        omms.print_energy_decomposition(simulation.context)
        if alch:
            log.info("Alchemical energy:")
            state = simulation.context.getState(getPositions=True)
            positions = state.getPositions()
            alch_ener = alch_factory.get_energy_components(alch_system, alch_state, positions)
            for key, val in alch_ener.items():
                print(key, val)

    # Calculate number of steps
    totstep = int(round(simultime/timestep, 0))
    log.info("Running for %s (%d steps)" % (simultime, totstep))
    log.info("Timestep is %s" % timestep)

    # Setup reporters
    log.info("Adding reporters...")
    ommlogfreq = int(round(outtime/timestep, 0))
    simulation.reporters.append(app.StateDataReporter(sys.stdout, ommlogfreq, step=True, time=True, 
            progress=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, 
            totalEnergy=True, temperature=True, speed=True, totalSteps=totstep, separator=' '))
    if ommlogfile:
        log.info("Adding OpenMM log file reporter to %s every %s (%d steps)" % (ommlogfile, outtime, ommlogfreq))
        simulation.reporters.append(app.StateDataReporter(ommlogfile, ommlogfreq, step=True, time=True, 
                progress=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, 
                totalEnergy=True, temperature=True, speed=True, totalSteps=totstep, separator=' '))
    if dcdoutfile:
        dcdoutfreq = int(round(outtime/timestep, 0))
        log.info("Adding DCD reporter to %s every %s (%d steps)" % (dcdoutfile, outtime, dcdoutfreq))
        simulation.reporters.append(app.DCDReporter(dcdoutfile, dcdoutfreq))
    if eneoutfile:
        eneoutfreq = int(round(outtime/timestep, 0))
        log.info("Adding Energy log file to %s every %s (%d steps)" % (eneoutfile, outtime, eneoutfreq))
        simulation.reporters.append(omms.EnergyReporter(eneoutfile, eneoutfreq))
    if alch and alch_outfile:
        alch_outfreq = int(round(outtime/timestep, 0))
        log.info("Adding Alchemical log file to %s every %s (%d steps)" % (alch_outfile, outtime, alch_outfreq))
        if alch_dynamics:
            simulation.reporters.append(omms.AlchemicalDynamics(alch_outfile, alch_outfreq, alch_system, 
                lambda_elec=alch_lambda_elec, lambda_vdw=alch_lambda_vdw, affect=alch_affect, 
                timestep=alch_dt, beta=alch_beta, dlambda=alch_dlambda))
        else:
            simulation.reporters.append(omms.AlchemicalReporter(alch_outfile, alch_outfreq, dlambda=alch_dlambda))
    if totstep>0 and not (dcdoutfile or ommlogfile or eneoutfile):
        log.warning("No reporter was added! You may want to set the <dcdoutfile>, <ommlogfile>, or <eneoutfile> variables!")

    #simulation.reporters.append(omms.CrashReporter("crash.pdb", savfreq))

    if totstep>0 and not xmloutfile:
        log.warning("The variable <xmloutfile> was not set. A restart file will not be written at the end of the simulation.")
    if totstep>0 and not pdboutfile:
        log.warning("The variable <pdboutfile> was not set. A pdb file with the final structure will not be written.")

    # Dynamics!
    if totstep<=0:
        log.info("Total number of steps negative or zero (%d). Do not run MD." % totstep)
    elif smt:
        log.info("Running simulated tempering")
        from simulated_tempering import SimulatedTempering
        st = SimulatedTempering(simulation, numTemperatures=15, 
                    minTemperature=300*u.kelvin, maxTemperature=450*u.kelvin, 
                    tempChangeInterval=25, reportInterval=ommlogfreq, 
                    reportFile="temp_"+dist+"_"+step+".dat")
        st.step(totstep)
    elif adt:
        log.info("Running adaptive tempering")
        from adaptive_tempering import AdaptiveTempering
        if restart==False:
            adtRstPrev=None
        AdaptiveTempering(simulation, totstep, temperature, mintemp=adtMinTemp, maxtemp=adtMaxTemp, timestep=1E-4, 
            stat_output_file=adtStatFile, output_file=adtRstFile, restart_file=adtRstPrev,
            inpfile=adtInpFile, logfile=adtLogFile,
            plugin_freq=adtPluginFreq, average_update_freq=adtAvgUpdt, temperature_update_freq=adtTempUpdt,
            stat_add_freq=ommlogfreq/adtPluginFreq, stat_output_freq=ommlogfreq/adtPluginFreq, output_file_freq=ommlogfreq/adtPluginFreq, 
            energy_interp_width=50, damping_average=adtDumping)
    else:
        log.info("Running plain dynamics for " + str(totstep) + " steps.")
        sys.stdout.flush()
        simulation.step(totstep)

    if xmloutfile:
        log.info("Saving status as restart file")
        simulation.saveState(xmloutfile)
    elif totstep>0:
        log.warning("No restart file was saved! Set <xmloutfile> variable next time!")

    #if pdboutfile:
    #    log.info("Writing pdb file")
    #    state = simulation.context.getState(getPositions=True)
    #    positions = state.getPositions()
    #    with open(pdboutfile, "w") as fp:
    #        pdb.writeModel(psf.topology, positions, fp)

    log.info("Normal termination :)")

    return simulation

