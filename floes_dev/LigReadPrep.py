#!/usr/bin/env python

# (C) 2021 OpenEye Scientific Software Inc. All rights reserved.
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

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.SubFloes.SubfloeFunctions import (setup_PLComplex_for_MD)

from MDOrion.Flask.cubes import ParallelRecordSizeCheck

from snowball import ExceptHandlerCube

job = WorkFloe('PL Complex formation for Short Trajectory MD',
               title='PL Complex formation for Short Trajectory MD')

job.description = """
PL Complex formation for the Short Trajectory MD (STMD) protocol based on
a prepared protein and a set of posed and prepared ligands as input.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
The protein needs to be prepared to MD standards: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
Given the inputs of the protein and posed ligands,
the complex is formed with each ligand/conformer
separately.

Required Input Parameters:
--------------------------
ligands (file): dataset of prepared ligands posed in the protein active site.
protein (file): dataset of the prepared protein structure.

Outputs:
--------
out (.oedb file): file of the protein-ligand complexes with parameters.
"""
# Locally the floe can be invoked by running the terminal command:
# python floes/LigReadPrep.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

floe_title = 'Protein-Ligand Flask Prep'
tags_for_floe = ['MDPrep']
#
job = WorkFloe(floe_title.join(' [{}]'.format(tag) for tag in tags_for_floe),
               title=floe_title)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = """
Starting with separate inputs for one MD-ready protein and multiple posed ligands,
this floe assembles a simulation-ready flask for each ligand, consisting of
the protein-ligand complex in a periodic box of explicit water with counterions.
If the protein input is not used, the floe will expect an MD-ready protein on the
same input record with each posed ligand, which it will use to form the complex with
that ligand.
"""

uuid = "cce33937-1eda-446a-864e-3627b58d09b4"

check_rec = ParallelRecordSizeCheck("Record Check Success")

exceptions = ExceptHandlerCube(floe_report_name="Analyze Floe Failure Report")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Out", description="MD Dataset out", order=2)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out", order=3)

job.add_cubes(check_rec, exceptions, ofs, fail)

# Call subfloe function to set up the solvated protein-ligand complex
PLComplex_for_MD_options = {}
PLComplex_for_MD_options['charge_ligands'] = True
PLComplex_for_MD_options['n_md_starts'] = 1
PLComplex_outcube = setup_PLComplex_for_MD(job, check_rec, PLComplex_for_MD_options)

# Connections
PLComplex_outcube.success.connect(check_rec.intake)
PLComplex_outcube.failure.connect(check_rec.fail_in)
check_rec.success.connect(ofs.intake)
check_rec.failure.connect(exceptions.intake)
exceptions.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
