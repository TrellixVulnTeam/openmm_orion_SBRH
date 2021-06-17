#!/usr/bin/env python

from floe.api import (WorkFloe, ParallelCubeGroup)

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ConformerGatheringData,
                                       ParallelConfTrajsToLigTraj,
                                       ParallelConcatenateTrajMMPBSACube)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                       ParallelMakeClusterTrajOEMols,
                                       ParallelMDTrajAnalysisClusterReport,
                                       ParallelClusterPopAnalysis,
                                       ParallelTrajAnalysisReportDataset,
                                       MDFloeReportCube)

job = WorkFloe("From combining confs Traj OEMols through Analysis")

job.description = """
Testing from aggregating conf trajs (into one ligand traj OEMol) through complete Analysis.
The input dataset is an .oedb file of the unaggregated confs MD results with Traj OEMols + IntE + PBSA
"""

# job.uuid = "7cacc2af-cae7-4dc7-8956-fcf539861e3d"

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="Flask Input OERecord", description="OERecord file name")

confGather = ConformerGatheringData("Gathering Conformer Records")
catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj")
catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube')
clusCube = ParallelClusterOETrajCube("ClusterOETrajCube")
clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis')
clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols')
prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset')
report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
report = MDFloeReportCube("report", title="Floe Report")

ofs = DatasetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="Flask Output OERecord", description="OERecord file name")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="Failures",
                       description="MD Dataset Failures out")

job.add_cubes(ifs,
              confGather, catLigTraj, catLigMMPBSA, clusCube, clusPop, clusOEMols,
              prepDataset, report_gen, report,
              ofs, fail)

ifs.success.connect(confGather.intake)
confGather.success.connect(catLigTraj.intake)
catLigTraj.success.connect(catLigMMPBSA.intake)
catLigMMPBSA.success.connect(clusCube.intake)
clusCube.success.connect(clusPop.intake)
clusPop.success.connect(clusOEMols.intake)
clusOEMols.success.connect(prepDataset.intake)
prepDataset.success.connect(report_gen.intake)
report_gen.success.connect(report.intake)
report_gen.failure.connect(fail.intake)
report.success.connect(ofs.intake)
report.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
