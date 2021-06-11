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

from MDOrion.SubFloes.SubfloeFunctions import setup_traj_analysis

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from snowball import (ExceptHandlerCube,
                      SuccessCounterCube)

from MDOrion.Flask.cubes import ParallelRecordSizeCheck

from MDOrion.Flask.cubes import CollectionSetting


floe_title = 'Analyze Protein-Ligand MD'
tags_for_floe = ['MDAnalysis']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = open(path.join(path.dirname(__file__), 'AnalyzePLMD_desc.rst'), 'r').read()
job.uuid = "7438db4d-30b1-478c-afc0-e921f0336c78"

# Ligand setting
iMDInput = DatasetReaderCube("MDInputReader", title="MD Input Reader")
iMDInput.promote_parameter("data_in", promoted_name="in",
                           title="MD Input Dataset", description="MD Input Dataset", order=0)

# This Cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenCollection", title="Open Collection")
coll_open.set_parameters(open=True)

# This Cube is necessary for the correct working of collection and shard
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

check_rec = ParallelRecordSizeCheck("Record Check Success", title="Record Check Success")

ofs = DatasetWriterCube('ofs', title='MD Analysis Output')
ofs.promote_parameter("data_out", promoted_name="out",
                      title="MD Analysis Output", description="MD Analysis output dataset name")

#ofs_du = DatasetWriterCube('ofs_du', title='MD Out Cluster Centric Output')
#ofs_du.promote_parameter("data_out", promoted_name="du",
#                         title="MD Out Cluster Centric Output", description="MD Out Cluster Centric Output")
exceptions = ExceptHandlerCube(floe_report_name="Analyze Floe Failure Report")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out", order=2)

job.add_cubes(iMDInput, coll_open, coll_close, check_rec,  ofs, exceptions, fail)

traj_anlys_outcube = setup_traj_analysis(job, coll_open, check_rec)

# Success Connections
iMDInput.success.connect(coll_open.intake)
traj_anlys_outcube.success.connect(coll_close.intake)
coll_close.success.connect(check_rec.intake)
check_rec.success.connect(ofs.intake)

# Fail Connections
traj_anlys_outcube.failure.connect(check_rec.fail_in)
coll_open.failure.connect(check_rec.fail_in)
coll_close.failure.connect(check_rec.fail_in)
check_rec.failure.connect(exceptions.intake)
exceptions.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
