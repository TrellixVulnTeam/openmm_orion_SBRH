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

from MDOrion.FEC.RFEC.cubes import (NESAnalysis,
                                    PredictDGFromDDG,
                                    PlotNESResults)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import MDFloeReportCube

from os import path

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.Flask.cubes import CollectionSetting

from MDOrion.Flask.cubes import ParallelRecordSizeCheck

from snowball import ExceptHandlerCube

floe_title = 'Non-Equilibrium Switching Recovery'
tags_for_floe = ['MD', 'FEC']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = open(path.join(path.dirname(__file__), 'NES_Recovery_desc.rst'), 'r').read()

job.uuid = "ebdc5d11-37e5-4f27-8752-a35170db03f9"

ines = DatasetReaderCube("NESRecoveryReader", title="NES Recovery Reader")
ines.promote_parameter("data_in", promoted_name="nes", title="NES Recovery Input Dataset", description="NES Recovery Input Dataset", order=0)

ibn = DatasetReaderCube("BoundReader", title="Bound Reader")
ibn.promote_parameter("data_in", promoted_name="bound", title="Bound Input Dataset", description="Bound Input Dataset", order=1)

nes_analysis = NESAnalysis("NES_Analysis")

ddg_to_dg = PredictDGFromDDG("RBFE to ABFE", title="RBFE to Affinity Estimate")
ddg_to_dg.promote_parameter('lig_exp_file', promoted_name='exp')

plot_aff = PlotNESResults("Plot Affinity Report", title="Plot Affinity Report")

report = MDFloeReportCube("NES Report", title="NES Floe Report")
report.set_parameters(floe_report_title="NES Report")

# This cube is necessary for the correct working of collections and shards
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

rec_check = ParallelRecordSizeCheck("Record Check Success", title="Record Size Checking")
rec_check_abfe = ParallelRecordSizeCheck("Record Check Success ABFE", title="Affinity Record Size Checking")

ofs_nes = DatasetWriterCube('ofs', title='NES Out')
ofs_nes.promote_parameter("data_out", promoted_name="out",
                          title="NES Dataset Out",
                          description="NES Dataset Out", order=2)

ofs_abfe = DatasetWriterCube('ofs_abfe', title='Affinity Out')
ofs_abfe.promote_parameter("data_out", promoted_name="abfe",
                           title="Affinity Out",
                           description="Affinity Out", order=3)

fail = DatasetWriterCube('fail', title='NES Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="NES Failures",
                       description="NES Dataset Failures out", order=4)

exceptions = ExceptHandlerCube(floe_report_name="Analyze Floe Failure Report")

job.add_cubes(ibn, ines, nes_analysis, ddg_to_dg, plot_aff, report, coll_close,
              rec_check_abfe, rec_check, exceptions, ofs_nes, ofs_abfe, fail)

ines.success.connect(nes_analysis.intake)
ibn.success.connect(ddg_to_dg.bound_port)

nes_analysis.success.connect(ddg_to_dg.intake)
nes_analysis.success.connect(report.intake)
ddg_to_dg.success.connect(rec_check_abfe.intake)
ddg_to_dg.graph_port.connect(plot_aff.intake)
plot_aff.success.connect(coll_close.intake)
coll_close.success.connect(rec_check.intake)

rec_check.success.connect(ofs_nes.intake)
rec_check_abfe.success.connect(ofs_abfe.intake)

# Fail Connections
nes_analysis.failure.connect(rec_check.fail_in)
ddg_to_dg.failure.connect(rec_check.fail_in)
report.failure.connect(rec_check.fail_in)
plot_aff.failure.connect(rec_check.fail_in)
coll_close.failure.connect(rec_check.fail_in)

rec_check.failure.connect(exceptions.intake)
rec_check_abfe.failure.connect(exceptions.intake)
exceptions.failure.connect(fail.intake)





