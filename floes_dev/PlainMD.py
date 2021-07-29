#!/usr/bin/env python

# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from os import path

from floe.api import (WorkFloe)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.Flask.cubes import MDComponentCube

from MDOrion.Flask.cubes import ParallelSolvationCube

from MDOrion.Flask.cubes import IDSettingCube

from MDOrion.ForceField.cubes import ForceFieldCube

from MDOrion.MDEngines.cubes import MDMinimizeCube

floe_title = 'Solvate and Run MD'
tags_for_floe = ['MDPrep', 'MDRun']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

ifs = DatasetReaderCube("SystemReader", title="Flask Reader")
ifs.promote_parameter("data_in", promoted_name="solute", title='Solute Input File',
                      description="Solute input file", order=0)

sysid = IDSettingCube("Flask Ids")

md_comp = MDComponentCube("MD Components")
md_comp.set_parameters(multiple_flasks=True)

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = ParallelSolvationCube("Hydration", title="Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')
solvate.set_parameters(close_solvent=True)

ff = ForceFieldCube("ForceField", title="Apply Force Field")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.1a1')

# Minimization
minComplex = MDMinimizeCube('minComplex', title='Minimization')
minComplex.modify_parameter(minComplex.restraints, promoted=False, default="noh (ligand or protein)")
minComplex.modify_parameter(minComplex.restraintWt, promoted=False, default=5.0)
minComplex.modify_parameter(minComplex.steps, promoted=False, default=0)
minComplex.set_parameters(center=True)
minComplex.set_parameters(save_md_stage=True)
minComplex.set_parameters(hmr=False)


ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out", order=1)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out", order=2)

job.add_cubes(ifs, sysid, md_comp, solvate, ff, minComplex, ofs, fail)

# Connections before setup_MD_startup subfloe
ifs.success.connect(sysid.intake)
sysid.success.connect(md_comp.intake)
md_comp.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
minComplex.success.connect(ofs.intake)

sysid.failure.connect(fail.intake)
md_comp.failure.connect(fail.intake)
solvate.failure.connect(fail.intake)
ff.failure.connect(fail.intake)
minComplex.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
