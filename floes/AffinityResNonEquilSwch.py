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

from floe.api import WorkFloe

from MDOrion.FEC.RFEC.cubes import PlotNESResults, PredictDGFromDDG

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from os import path

from MDOrion.Flask.cubes import ParallelRecordSizeCheck


floe_title = 'Compare Experimental Affinity with NES Results'
tags_for_floe = ['Utility', 'FECalc']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = open(path.join(path.dirname(__file__), 'AffinityResNonEquilSwch_desc.rst'), 'r').read()

job.uuid = "e50c2e49-63a0-4bb3-aba7-b65dc4f92ab9"

ifs = DatasetReaderCube("NESInputReader", title="NES Input Reader")
ifs.promote_parameter("data_in", promoted_name="plot",
                      title='NES Input Reader',
                      description="The Dataset produced by the Non-Equilibrium Switching", order=0)

bnd_eq = DatasetReaderCube("BoundEqReader", title="Bound Equilibrium Reader")
bnd_eq.promote_parameter("data_in", promoted_name="bound",
                         title='Bound Equilibrium Reader',
                         description="The Equilibrium Bound Dataset")

ddg_to_dg_sub = PredictDGFromDDG("RBFE to Affinity Estimate", title="RBFE to Affinity Estimate")
ddg_to_dg_sub.promote_parameter('lig_exp_file', promoted_name='exp', required=True)

plot = PlotNESResults("PlotAffinities", title="Plot Affinities")

ofs_DG = DatasetWriterCube('ofs_DG', title='Affinity Out')
ofs_DG.promote_parameter("data_out", promoted_name="DG",
                         title="Affinity Out",
                         description="Affinity Out", order=1)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail",
                       title="Failures",
                       description="Fail Data Set", order=2)


rec_check = ParallelRecordSizeCheck("Record Check Success", title="Record Size Checking")


job.add_cubes(ifs, bnd_eq, ddg_to_dg_sub, plot, fail, rec_check, ofs_DG)

ifs.success.connect(ddg_to_dg_sub.intake)

bnd_eq.success.connect(ddg_to_dg_sub.bound_port)
ddg_to_dg_sub.graph_port.connect(plot.intake)
ddg_to_dg_sub.success.connect(rec_check.intake)
rec_check.success.connect(ofs_DG.intake)

ddg_to_dg_sub.failure.connect(rec_check.fail_in)
plot.failure.connect(rec_check.fail_in)

rec_check.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
