#############################################################################
# Copyright (C) 2018 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np
from openeye import oechem
from openeye import oezap
from simtk import (unit, openmm)

def ParmedAndOEMolAtomsMatch( parmedObj, mol):
    # first verify that the total number of atoms match
    pmdAtoms = parmedObj.atoms
    if len(pmdAtoms)!=mol.NumAtoms():
        print('Unequal number of atoms between parmed ({}) and OEMol({})'.
             format(len(pmdAtoms),mol.NumAtoms()))
        return False
    # second verify each atom's atomic numbers match
    for pmdAtm, oeAtom in zip(pmdAtoms,mol.GetAtoms()):
        if pmdAtm.atomic_number!=oeAtom.GetAtomicNum():
            print('different atomic number between parmed ({}) and OEMol({}) atom {}'.
                  format(pmdAtm.atomic_number,oeAtom.GetAtomicNum(),pmdAtm.name))
            return False
    #
    return True


def ProtLigInteractionEFromParmedOETraj( pmed, ligOETraj, protOETraj):
    # Generate ligand and protein parmed subsystems
    print('Generating ligand and protein parmed subsystems')
    proteinPmed = pmed.split()[0][0]
    ligandPmed = pmed.split()[1][0]
    complexPmed = proteinPmed + ligandPmed
    #complexPmed.save("system.pdb", overwrite=True)
    #
    # Check that protein and ligand atoms match between parmed subsystems, ligOETraj, and protOETraj
    if not ParmedAndOEMolAtomsMatch( ligandPmed, ligOETraj):
        print('ligand atoms do not match between parmed and ligOETraj')
        return None, None, None, None
    if not ParmedAndOEMolAtomsMatch( proteinPmed, protOETraj):
        print('protein atoms do not match between parmed and protOETraj')
        return None, None, None, None
    print('protein and ligand atoms match between parmed subsystems, ligOETraj, and protOETraj')
    #
    # While we are here put the parmed charges on protOETraj and ligOETraj as side effects
    for pmdAtom, oeAtom in zip( proteinPmed.atoms, protOETraj.GetAtoms()):
        oeAtom.SetPartialCharge( pmdAtom.charge)
    for pmdAtom, oeAtom in zip( ligandPmed.atoms, ligOETraj.GetAtoms()):
        oeAtom.SetPartialCharge( pmdAtom.charge)
    #
    # Generate openmm components needed for energy evals
    # ligand
    ligOmmSys = ligandPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
    ligIntegrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
                                              0.002 * unit.picoseconds)
    ligSim = openmm.app.Simulation(ligandPmed.topology, ligOmmSys, ligIntegrator)
    # protein
    protOmmSys = proteinPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
    protIntegrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
                                              0.002 * unit.picoseconds)
    protSim = openmm.app.Simulation(proteinPmed.topology, protOmmSys, protIntegrator)
    # complex
    cplxOmmSys = complexPmed.createSystem(nonbondedMethod=openmm.app.NoCutoff)
    cplxIntegrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
                                              0.002 * unit.picoseconds)
    cplxSim = openmm.app.Simulation(complexPmed.topology, cplxOmmSys, cplxIntegrator)
    print('OpenMM components generated for protein, ligand, and complex subsystems')
    #
    # Evaluate energies for the subsystems
    ligE = []
    protE = []
    cplxE = []
    plIntE = []
    #
    for ligConf, protConf in zip(ligOETraj.GetConfs(), protOETraj.GetConfs()):
        ligXYZdict = ligConf.GetCoords()
        ligXYZ = [openmm.Vec3(v[0], v[1], v[2]) for k, v in ligXYZdict.items()] * unit.angstrom
        ligSim.context.setPositions(ligXYZ)
        ligState = ligSim.context.getState(getEnergy=True)
        ligIntraE = ligState.getPotentialEnergy().in_units_of(
            unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
        ligE.append( ligIntraE)

        protXYZdict = protConf.GetCoords()
        protXYZ = [openmm.Vec3(v[0], v[1], v[2]) for k, v in protXYZdict.items()] * unit.angstrom
        protSim.context.setPositions(protXYZ)
        protState = protSim.context.getState(getEnergy=True)
        protIntraE = protState.getPotentialEnergy().in_units_of(
            unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
        protE.append( protIntraE)

        cplxXYZ = protXYZ + ligXYZ
        cplxSim.context.setPositions(cplxXYZ)
        cplxState = cplxSim.context.getState(getEnergy=True)
        cplxIntraE = cplxState.getPotentialEnergy().in_units_of(
            unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
        cplxE.append( cplxIntraE)

        plIntE.append( cplxIntraE-(protIntraE+ligIntraE) )

    print('OpenMM energies computed for protein, ligand, and complex trajectories')
    return plIntE, cplxE, protE, ligE


def PBSA( ligand, protein):
    bind = oezap.OEBind()
    bind.SetProtein(protein)
    results = oezap.OEBindResults()
    if not bind.Bind(ligand, results):
        print( 'zap Bind run failed for {} with {}'.format( ligand.GetTitle(), protein.GetTitle()))
        return None, None, None, None, None, None
    # convert key values from kT to kcal/mol
    kTtoKcal = 0.59
    Ebind = kTtoKcal*results.GetBindingEnergy()
    EbindEl = kTtoKcal*results.GetZapEnergy()
    EdesolEl = kTtoKcal*results.GetDesolvationEnergy()
    EintEl = kTtoKcal*( results.GetZapEnergy()-results.GetDesolvationEnergy() )
    EbindSA = kTtoKcal*results.GetBuriedAreaEnergy()
    SAburied = results.GetBuriedArea()
    #
    return Ebind, EbindEl, EdesolEl, EintEl, EbindSA, SAburied


def TrajPBSA( ligOETraj, protOETraj, radiiType=oechem.OERadiiType_Zap9):
    # set ZAP9 radii on protein and ligand
    oechem.OEAssignRadii(protOETraj, radiiType, oechem.OERadiiType_HonigIonicCavity)
    oechem.OEAssignRadii(ligOETraj, radiiType)
    #
    zapBind = []
    zapBindEl = []
    zapDesolEl = []
    zapIntEl = []
    zapBindSA25 = []
    saBuried = []
    #
    for protConf, ligConf in zip( protOETraj.GetConfs(), ligOETraj.GetConfs() ):
        Ebind, EbindEl, EdesolEl, EintEl, EbindSA, buriedSA = PBSA( ligConf, protConf)
        zapBind.append( Ebind)
        zapBindEl.append( EbindEl)
        zapDesolEl.append( EdesolEl)
        zapIntEl.append( EintEl)
        zapBindSA25.append( EbindSA)
        saBuried.append( buriedSA)
    #
    return zapBind, zapBindEl, zapDesolEl, zapIntEl, zapBindSA25, saBuried
