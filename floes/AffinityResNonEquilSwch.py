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

job = WorkFloe("Plot the Affinity results for NES",
               title="Plot the Affinity results for the Non-Equilibrium Switching floe")


job.description = open(path.join(path.dirname(__file__), 'AffinityResNonEquilSwch.rst'), 'r').read()

job.classification = [['NES Results']]
job.uuid = "e50c2e49-63a0-4bb3-aba7-b65dc4f92ab9"
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="Plot Reader")
ifs.promote_parameter("data_in", promoted_name="plot",
                      title='Plot Input File',
                      description="The Dataset produced by the Short Trajectory MD with Analysis floe", order=0)

ddg_to_dg_sub = PredictDGFromDDG("RBFE to ABFE", title="RBFE to Affinity Estimate")
ddg_to_dg_sub.promote_parameter('lig_exp_file', promoted_name='exp', required=True)

plot = PlotNESResults("PlotAffinities", title="Plot Affinities")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", description="Fail Data Set", order=2)

rec_check = ParallelRecordSizeCheck("Record Check Success", title="Record Size Checking")

job.add_cubes(ifs, ddg_to_dg_sub, plot, fail, rec_check)

ifs.success.connect(ddg_to_dg_sub.intake)
ddg_to_dg_sub.graph_port.connect(plot.intake)

ddg_to_dg_sub.failure.connect(rec_check.fail_in)
plot.failure.connect(rec_check.fail_in)

rec_check.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
