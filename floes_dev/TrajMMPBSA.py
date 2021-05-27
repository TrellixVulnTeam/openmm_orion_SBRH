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

from MDOrion.SubFloes.SubfloeFunctions import (setup_traj_mmpbsa)

from MDOrion.Flask.cubes import ParallelRecordSizeCheck

from snowball import ExceptHandlerCube


floe_title = 'MMPBSA Processing of Protein-Ligand MD Dataset'
tags_for_floe = ['MDAnalysis']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = """
An MMPBSA analysis is carried out on trajectories from a Protein-Ligand MD Dataset.

Required Input Parameters:
--------------------------
.oedb of records from a Protein-Ligand MD Dataset.

Outputs:
--------
out (.oedb file): Dataset of the Traj OEMols and MMPBSA  results.
"""

job.uuid = "2717cf39-5bdd-4a1e-880e-5208bb232959"


ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="Flask Input OERecord", description="OERecord file name")

check_rec = ParallelRecordSizeCheck("Record Check Success")

exceptions = ExceptHandlerCube(floe_report_name="Analyze Floe Failure Report")

ofs = DatasetWriterCube('ofs', title='MD Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs, check_rec, exceptions, ofs, fail)

# Call subfloe function to process the trajectories
traj_proc_outcube = setup_traj_mmpbsa(job, ifs, check_rec)

# Connections
traj_proc_outcube.success.connect(check_rec.intake)
traj_proc_outcube.failure.connect(check_rec.fail_in)
check_rec.success.connect(ofs.intake)
check_rec.failure.connect(exceptions.intake)
exceptions.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
