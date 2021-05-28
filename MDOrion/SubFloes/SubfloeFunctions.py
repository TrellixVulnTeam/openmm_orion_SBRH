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

from floe.api import ParallelCubeGroup

from orionplatform.cubes import DatasetReaderCube


from MDOrion.Flask.cubes import (IDSettingCube,
                                 ParallelSolvationCube,
                                 MDComponentCube)

from MDOrion.MDEngines.cubes import (ParallelMDMinimizeCube,
                                     ParallelMDNvtCube,
                                     ParallelMDNptCube)

from MDOrion.ComplexPrep.cubes import ComplexPrepCube

from MDOrion.ForceField.cubes import ParallelForceFieldCube

from MDOrion.LigPrep.cubes import (ParallelLigandChargeCube,
                                   LigandSetting)

from MDOrion.FEC.RFEC.cubes import (BoundUnboundSwitchCube,
                                    RBFECEdgeGathering,
                                    ParallelNESGMXChimera,
                                    ParallelNESGMX,
                                    NESAnalysis,
                                    PredictDGFromDDG,
                                    PlotNESResults)

from MDOrion.TrjAnalysis.cubes_trajProcessing import (ParallelTrajToOEMolCube,
                                                      ParallelTrajInteractionEnergyCube,
                                                      ParallelTrajPBSACube,
                                                      ConformerGatheringData,
                                                      ParallelConfTrajsToLigTraj,
                                                      ParallelConcatenateTrajMMPBSACube)

from MDOrion.TrjAnalysis.cubes_clusterAnalysis import (ParallelClusterOETrajCube,
                                                       ParallelMakeClusterTrajOEMols,
                                                       ParallelMDTrajAnalysisClusterReport,
                                                       ParallelClusterPopAnalysis,
                                                       ParallelTrajAnalysisReportDataset,
                                                       MDFloeReportCube,)

from MDOrion.TrjAnalysis.cubes_hintAnalysis import (ParallelBintScoreInitialPoseAndTrajectory)

from snowball.utils.dataset_reader_opt import DatasetReaderOptCube


def nes_gmx_subfloe(floe_job, input_port_dic, output_port_dic, options):

    input_bound_port = input_port_dic['input_bound_port']
    input_open_collection_port = input_port_dic['input_open_collection_port']

    output_nes_port = output_port_dic['output_nes_port']
    output_DG_port = output_port_dic['output_DG_port']
    output_recovery_port = output_port_dic['output_recovery']
    output_fail_port = output_port_dic['output_fail_port']

    # Switching Bound and Unbound runs
    switch_sub = BoundUnboundSwitchCube("Bound/Unbound Switch NES", title='Bound/Unbound Switch NES')

    gathering_sub = RBFECEdgeGathering("Gathering", title="Gathering Equilibrium Runs")
    gathering_sub.promote_parameter('map_file', promoted_name=options['edge_map_file'], order=2)

    chimera_sub = ParallelNESGMXChimera("GMXChimera", title="GMX Chimera")
    chimera_sub.promote_parameter("trajectory_frames", promoted_name="trajectory_frames",
                                  default=options['n_traj_frames'],
                                  description="The total number of trajectory frames to be used along the NE switching", order=2)

    unbound_nes_sub = ParallelNESGMX("GMXUnboundNES", title="GMX Unbound NES")
    unbound_nes_sub.promote_parameter("time", promoted_name="nes_time",
                                      default=options['nes_switch_time_in_ns'], order=3)

    unbound_nes_sub.modify_parameter(unbound_nes_sub.instance_type, promoted=False, default='c5')
    unbound_nes_sub.modify_parameter(unbound_nes_sub.cpu_count, promoted=False, default=2)
    unbound_nes_sub.modify_parameter(unbound_nes_sub.gpu_count, promoted=False, default=0)
    unbound_nes_sub.modify_parameter(unbound_nes_sub.memory_mb, promoted=False, default=3*1024)
    unbound_nes_sub.modify_parameter(unbound_nes_sub.gmx_mpi_threads, promoted=False, default=1)
    unbound_nes_sub.modify_parameter(unbound_nes_sub.gmx_openmp_threads, promoted=False, default=2)
    unbound_nes_sub.modify_parameter(unbound_nes_sub.max_parallel, promoted=False, default=10000)

    bound_nes_sub = ParallelNESGMX("GMXBoundNES", title="GMX Bound NES")
    bound_nes_sub.promote_parameter("time", promoted_name="nes_time", order=3)
    # bound_nes.modify_parameter(bound_nes.instance_type, promoted=False, default='g3')
    bound_nes_sub.modify_parameter(bound_nes_sub.cpu_count, promoted=False, default=16)
    bound_nes_sub.modify_parameter(bound_nes_sub.gpu_count, promoted=False, default=1)
    bound_nes_sub.modify_parameter(bound_nes_sub.gmx_mpi_threads, promoted=False, default=1)
    bound_nes_sub.modify_parameter(bound_nes_sub.gmx_openmp_threads, promoted=False, default=16)
    bound_nes_sub.modify_parameter(bound_nes_sub.max_parallel, promoted=False, default=10000)

    nes_analysis_sub = NESAnalysis("NES_Analysis")

    ddg_to_dg_sub = PredictDGFromDDG("RBFE to ABFE", title="RBFE to Affinity Estimate")
    ddg_to_dg_sub.promote_parameter('lig_exp_file', promoted_name='exp')

    plot_aff_sub = PlotNESResults("Plot Affinity Report", title="Plot Affinity Report")

    report_sub = MDFloeReportCube("NES Report", title="NES Floe Report")
    report_sub.set_parameters(floe_report_title="NES Report")

    floe_job.add_cubes(switch_sub, gathering_sub,
                       chimera_sub, bound_nes_sub, unbound_nes_sub,
                       nes_analysis_sub, ddg_to_dg_sub,
                       plot_aff_sub, report_sub)

    input_open_collection_port.connect(switch_sub.intake)

    switch_sub.success.connect(gathering_sub.intake)
    switch_sub.bound_port.connect(gathering_sub.bound_port)

    # Chimera NES Setting
    gathering_sub.success.connect(chimera_sub.intake)

    chimera_sub.bound_port.connect(bound_nes_sub.intake)
    bound_nes_sub.success.connect(nes_analysis_sub.intake)

    chimera_sub.success.connect(unbound_nes_sub.intake)
    unbound_nes_sub.success.connect(nes_analysis_sub.intake)

    bound_nes_sub.success.connect(output_recovery_port)
    unbound_nes_sub.success.connect(output_recovery_port)

    nes_analysis_sub.success.connect(ddg_to_dg_sub.intake)
    input_bound_port.connect(ddg_to_dg_sub.bound_port)

    ddg_to_dg_sub.success.connect(output_DG_port)

    ddg_to_dg_sub.graph_port.connect(report_sub.intake)
    report_sub.success.connect(plot_aff_sub.intake)
    plot_aff_sub.success.connect(output_nes_port)

    # Fail port connections
    switch_sub.failure.connect(output_fail_port)

    gathering_sub.failure.connect(output_fail_port)
    chimera_sub.failure.connect(output_fail_port)
    unbound_nes_sub.failure.connect(output_fail_port)
    bound_nes_sub.failure.connect(output_fail_port)

    nes_analysis_sub.failure.connect(output_fail_port)
    ddg_to_dg_sub.failure.connect(output_fail_port)
    plot_aff_sub.failure.connect(output_fail_port)
    report_sub.failure.connect(output_fail_port)

    return True


def setup_PLComplex_for_MD(input_floe, fail_cube, options):

    # Ligand setting
    iligs = DatasetReaderCube("LigandReader", title="Ligand Reader")
    iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input Dataset",
                            description="Ligand Dataset", order=1)

    ligset = LigandSetting("LigandSetting", title="Ligand Setting")
    ligset.promote_parameter('max_md_runs', promoted_name='max_md_runs',
                             default=500,
                             description='The maximum allowed number of md runs')
    ligset.promote_parameter('n_md_starts', promoted_name='n_md_starts',
                             default=options['n_md_starts'],
                             description='The number of md starts for each ligand/conformer')
    ligset.set_parameters(lig_res_name='LIG')

    chargelig = ParallelLigandChargeCube("LigCharge", title="Ligand Charge")
    chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                                description="Charge the ligand or not", default=options['charge_ligands'])

    ligid = IDSettingCube("Ligand Ids")

    # Protein Reading cube. The protein prefix parameter is used to select a name for the
    # output system files
    iprot = DatasetReaderOptCube("ProteinReader", title="Protein Reader")
    iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input Dataset',
                            description="Protein Dataset", order=0)

    # Protein Setting
    mdcomp = MDComponentCube("MD Components", title="MD Components")
    mdcomp.promote_parameter("flask_title", promoted_name="flask_title", default="")

    # Complex cube used to assemble the ligands and the solvated protein
    complx = ComplexPrepCube("Complex", title="Complex Preparation")

    # The solvation cube is used to solvate the system and define the ionic strength of the solution
    solvate = ParallelSolvationCube("Solvation", title="Solvation")

    input_floe.add_cubes(iligs, ligset, chargelig, ligid, iprot, mdcomp, complx, solvate)

    # Success Connections
    iligs.success.connect(ligset.intake)
    ligset.success.connect(chargelig.intake)
    chargelig.success.connect(ligid.intake)
    ligid.success.connect(complx.intake)
    iprot.success.connect(mdcomp.intake)
    mdcomp.success.connect(complx.protein_port)
    complx.success.connect(solvate.intake)

    # Fail Connections
    ligset.failure.connect(fail_cube.fail_in)
    chargelig.failure.connect(fail_cube.fail_in)
    ligid.failure.connect(fail_cube.fail_in)
    mdcomp.failure.connect(fail_cube.fail_in)
    complx.failure.connect(fail_cube.fail_in)

    return solvate


def setup_MD_startup(input_floe, input_cube, fail_cube, options):
    # Force Field Application
    ff = ParallelForceFieldCube("ForceField", title="Apply Force Field")
    ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='Amber14SB')
    ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.0')

    # Production run
    prod = ParallelMDNptCube("Production", title="Production")
    prod.promote_parameter('time', promoted_name='prod_ns',
                           default=options['Prod_Default_Time_ns'],
                           description='Length of MD run in nanoseconds')
    prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval',
                           default=options['Prod_Default_Traj_Intvl_ns'],
                           description='Trajectory saving interval in ns')
    prod.promote_parameter('hmr', promoted_name="HMR", title='Use Hydrogen Mass Repartitioning', default=True,
                           description='Give hydrogens more mass to speed up the MD')
    prod.promote_parameter('md_engine', promoted_name='md_engine', default='OpenMM',
                           description='Select the MD Engine')
    prod.set_parameters(reporter_interval=options['Prod_Default_Traj_Intvl_ns'])
    prod.set_parameters(suffix='prod')

    # Minimization
    minComplex = ParallelMDMinimizeCube('minComplex', title='Minimization')
    minComplex.modify_parameter(minComplex.restraints, promoted=False, default="noh (ligand or protein)")
    minComplex.modify_parameter(minComplex.restraintWt, promoted=False, default=5.0)
    minComplex.modify_parameter(minComplex.steps, promoted=False, default=0)
    minComplex.set_parameters(center=True)
    minComplex.set_parameters(save_md_stage=True)
    minComplex.set_parameters(hmr=False)
    minComplex.promote_parameter("md_engine", promoted_name="md_engine")

    # NVT simulation. Here the assembled system is warmed up to the final selected temperature
    warmup = ParallelMDNvtCube('warmup', title='Warm Up')
    warmup.set_parameters(time=0.01)
    warmup.modify_parameter(warmup.restraints, promoted=False, default="noh (ligand or protein)")
    warmup.modify_parameter(warmup.restraintWt, promoted=False, default=2.0)
    warmup.set_parameters(trajectory_interval=0.0)
    warmup.set_parameters(reporter_interval=0.001)
    warmup.set_parameters(suffix='warmup')
    warmup.set_parameters(hmr=False)
    warmup.set_parameters(save_md_stage=True)
    warmup.promote_parameter("md_engine", promoted_name="md_engine")

    # The system is equilibrated at the right pressure and temperature in several stages
    # The main difference between the stages is related to the restraint force used
    # to keep the ligand and protein in their starting positions. A relatively strong force
    # is applied in the first stage while a relatively small one is applied in the latter

    # NPT Equilibration stage 1
    equil1 = ParallelMDNptCube('equil1', title='Equilibration I')
    equil1.set_parameters(time=0.01)
    equil1.promote_parameter("hmr", promoted_name="HMR", default=True)
    equil1.modify_parameter(equil1.restraints, promoted=False, default="noh (ligand or protein)")
    equil1.modify_parameter(equil1.restraintWt, promoted=False, default=1.0)
    equil1.set_parameters(trajectory_interval=0.0)
    equil1.set_parameters(reporter_interval=0.001)
    equil1.set_parameters(suffix='equil1')
    equil1.promote_parameter("md_engine", promoted_name="md_engine")

    # NPT Equilibration stage 2
    equil2 = ParallelMDNptCube('equil2', title='Equilibration II')
    equil2.set_parameters(time=0.02)
    equil2.promote_parameter("hmr", promoted_name="HMR", default=True)
    equil2.modify_parameter(equil2.restraints, promoted=False, default="noh (ligand or protein)")
    equil2.modify_parameter(equil2.restraintWt, promoted=False, default=0.5)
    equil2.set_parameters(trajectory_interval=0.0)
    equil2.set_parameters(reporter_interval=0.001)
    equil2.set_parameters(suffix='equil2')
    equil2.promote_parameter("md_engine", promoted_name="md_engine")

    # NPT Equilibration stage 3
    equil3 = ParallelMDNptCube('equil3', title='Equilibration III')
    equil3.modify_parameter(equil3.time, promoted=False, default=0.1)
    equil3.promote_parameter("hmr", promoted_name="HMR", default=True)
    equil3.modify_parameter(equil3.restraints, promoted=False, default="noh (ligand or protein)")
    equil3.modify_parameter(equil3.restraintWt, promoted=False, default=0.2)
    equil3.set_parameters(trajectory_interval=0.0)
    equil3.set_parameters(reporter_interval=0.002)
    equil3.set_parameters(suffix='equil3')
    equil3.promote_parameter("md_engine", promoted_name="md_engine")

    # NPT Equilibration stage 4
    equil4 = ParallelMDNptCube('equil4', title='Equilibration IV')
    equil4.modify_parameter(equil4.time, promoted=False, default=0.1)
    equil4.promote_parameter("hmr", promoted_name="HMR", default=True)
    equil4.modify_parameter(equil4.restraints, promoted=False, default="ca_protein or (noh ligand)")
    equil4.modify_parameter(equil4.restraintWt, promoted=False, default=0.1)
    equil4.set_parameters(trajectory_interval=0.0)
    equil4.set_parameters(reporter_interval=0.002)
    equil4.set_parameters(suffix='equil4')
    equil4.promote_parameter("md_engine", promoted_name="md_engine")

    md_group = ParallelCubeGroup(cubes=[minComplex, warmup, equil1, equil2, equil3, equil4, prod])
    input_floe.add_group(md_group)

    input_floe.add_cubes(ff, minComplex, warmup, equil1, equil2, equil3, equil4, prod)

    # Success Connections
    input_cube.success.connect(ff.intake)
    ff.success.connect(minComplex.intake)
    minComplex.success.connect(warmup.intake)
    warmup.success.connect(equil1.intake)
    equil1.success.connect(equil2.intake)
    equil2.success.connect(equil3.intake)
    equil3.success.connect(equil4.intake)
    equil4.success.connect(prod.intake)
    
    # Fail Connections
    ff.failure.connect(fail_cube.fail_in)
    minComplex.failure.connect(fail_cube.fail_in)
    warmup.failure.connect(fail_cube.fail_in)
    equil1.failure.connect(fail_cube.fail_in)
    equil2.failure.connect(fail_cube.fail_in)
    equil3.failure.connect(fail_cube.fail_in)
    equil4.failure.connect(fail_cube.fail_in)

    return prod


def setup_MDsmallmol_startup(input_floe, input_cube, fail_cube, options):
    # Force Field Application
    ff_small = ParallelForceFieldCube("ForceField", title="Apply Force Field")
    ff_small.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='OpenFF_1.3.0')

    # Production run
    prod_small = ParallelMDNptCube("Production Unbound State", title="Production Unbound State")
    prod_small.promote_parameter('time', promoted_name='prod_unb_ns',
                                 default=options['Prod_Default_Time_ns'],
                                 description='Length of Unbound MD run in nanoseconds')
    prod_small.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval',
                           default=options['Prod_Default_Traj_Intvl_ns'],
                           description='Unbound trajectory saving interval in ns')
    prod_small.promote_parameter('hmr', promoted_name="hmr_us",
                                 title='Use Hydrogen Mass Repartitioning in the Unbound simulation',
                                 default=True,
                                 description='Give hydrogens more mass to speed up the MD')
    prod_small.set_parameters(md_engine='OpenMM')
    prod_small.set_parameters(reporter_interval=options['Prod_Default_Traj_Intvl_ns'])
    prod_small.set_parameters(suffix='prod_unb')

    # Minimization of the Unbound-States
    minimize_small = ParallelMDMinimizeCube("Minimize Unbound States", title="Minimization Unbound States")
    minimize_small.set_parameters(restraints='noh ligand')
    minimize_small.set_parameters(md_engine='OpenMM')
    minimize_small.set_parameters(steps=2000)
    minimize_small.set_parameters(restraintWt=5.0)
    minimize_small.set_parameters(center=True)
    minimize_small.set_parameters(hmr=False)

    # NVT Warm-up of the Unbound-States
    warmup_small = ParallelMDNvtCube('Warmup Unbound States', title='Warmup Unbound States')
    warmup_small.set_parameters(time=0.01)
    warmup_small.set_parameters(restraints="noh ligand")
    warmup_small.set_parameters(md_engine='OpenMM')
    warmup_small.set_parameters(restraintWt=2.0)
    warmup_small.set_parameters(trajectory_interval=0.0)
    warmup_small.set_parameters(reporter_interval=0.002)
    warmup_small.set_parameters(suffix='warmup_unb')
    warmup_small.set_parameters(hmr=False)

    # NPT Equilibration stage of the Unbound-States
    equil_small = ParallelMDNptCube('Equilibration Unbound State', title='Equilibration Unbound States')
    equil_small.set_parameters(time=0.1)
    equil_small.promote_parameter("hmr", promoted_name="hmr_unbound", default=True)
    equil_small.set_parameters(restraints="noh ligand")
    equil_small.set_parameters(md_engine='OpenMM')
    equil_small.set_parameters(restraintWt=0.1)
    equil_small.set_parameters(trajectory_interval=0.0)
    equil_small.set_parameters(reporter_interval=0.004)
    equil_small.set_parameters(suffix='equil_unb')

    md_group_small = ParallelCubeGroup(cubes=[minimize_small, warmup_small, equil_small, prod_small])
    input_floe.add_group(md_group_small)

    input_floe.add_cubes(ff_small, minimize_small, warmup_small, equil_small, prod_small)

    # Success port Connections
    input_cube.success.connect(ff_small.intake)
    ff_small.success.connect(minimize_small.intake)
    minimize_small.success.connect(warmup_small.intake)
    warmup_small.success.connect(equil_small.intake)
    equil_small.success.connect(prod_small.intake)

    # Fail port connections
    input_cube.failure.connect(fail_cube.intake)
    ff_small.failure.connect(fail_cube.intake)
    minimize_small.failure.connect(fail_cube.intake)
    warmup_small.failure.connect(fail_cube.intake)
    equil_small.failure.connect(fail_cube.intake)

    return prod_small


def setup_traj_mmpbsa(input_floe, input_cube, fail_cube):
    trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube", title="Trajectory To OEMols")
    IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube", title="MM Energies")
    PBSACube = ParallelTrajPBSACube("TrajPBSACube", title="PBSA Energies")

    trajproc_group = ParallelCubeGroup(cubes=[trajCube, IntECube, PBSACube])
    input_floe.add_group(trajproc_group)

    input_floe.add_cubes(trajCube, IntECube, PBSACube)

    # Success Connections
    input_cube.success.connect(trajCube.intake)
    trajCube.success.connect(IntECube.intake)
    IntECube.success.connect(PBSACube.intake)

    # Fail Connections
    trajCube.failure.connect(fail_cube.fail_in)
    IntECube.failure.connect(fail_cube.fail_in)

    return PBSACube


def setup_gather_cluster(input_floe, input_cube, fail_cube):
    confGather = ConformerGatheringData("Gathering Conformer Records", title="Gathering Conformer Records")
    catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj", title="Combine Pose Trajectories")
    catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube', title="Concatenate MMPBSA Energies")
    trajBints = ParallelBintScoreInitialPoseAndTrajectory("TrajBintsCube", title="Trajectory Binding Interactions")
    clusCube = ParallelClusterOETrajCube("ClusterOETrajCube", title="Clustering")
    clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis', title="Clustering Analysis")
    clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols', title="Per-Cluster Analysis")
    prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset', title="Analysis Report")
    report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport", title="Relevant Output Extraction")

    analysis_group = ParallelCubeGroup(cubes=[catLigTraj, catLigMMPBSA, trajBints, clusCube, clusPop,
                                              clusOEMols, prepDataset, report_gen])
    input_floe.add_group(analysis_group)

    report = MDFloeReportCube("report", title="Floe Report")

    input_floe.add_cubes(confGather,
                  catLigTraj, catLigMMPBSA, trajBints, clusCube, clusPop, clusOEMols,
                  prepDataset, report_gen, report)

    # Success Connections
    input_cube.success.connect(confGather.intake)
    confGather.success.connect(catLigTraj.intake)
    catLigTraj.success.connect(catLigMMPBSA.intake)
    catLigMMPBSA.success.connect(trajBints.intake)
    trajBints.success.connect(clusCube.intake)
    clusCube.success.connect(clusPop.intake)
    clusPop.success.connect(clusOEMols.intake)
    clusOEMols.success.connect(prepDataset.intake)
    prepDataset.success.connect(report_gen.intake)
    report_gen.success.connect(report.intake)

    # Fail Connections
    confGather.failure.connect(fail_cube.fail_in)
    catLigTraj.failure.connect(fail_cube.fail_in)
    catLigMMPBSA.failure.connect(fail_cube.fail_in)
    trajBints.failure.connect(fail_cube.fail_in)
    clusCube.failure.connect(fail_cube.fail_in)
    clusPop.failure.connect(fail_cube.fail_in)
    clusOEMols.failure.connect(fail_cube.fail_in)
    prepDataset.failure.connect(fail_cube.fail_in)
    report_gen.failure.connect(fail_cube.fail_in)

    return report



def setup_traj_analysis(input_floe, input_cube, fail_cube, du_output_cube):
    trajCube = ParallelTrajToOEMolCube("TrajToOEMolCube", title="Trajectory To OEMols")
    IntECube = ParallelTrajInteractionEnergyCube("TrajInteractionEnergyCube", title="MM Energies")
    PBSACube = ParallelTrajPBSACube("TrajPBSACube", title="PBSA Energies")

    trajproc_group = ParallelCubeGroup(cubes=[trajCube, IntECube, PBSACube])
    input_floe.add_group(trajproc_group)

    confGather = ConformerGatheringData("Gathering Conformer Records", title="Gathering Conformer Records")
    catLigTraj = ParallelConfTrajsToLigTraj("ConfTrajsToLigTraj", title="Combine Pose Trajectories")
    catLigMMPBSA = ParallelConcatenateTrajMMPBSACube('ConcatenateTrajMMPBSACube', title="Concatenate MMPBSA Energies")
    trajBints = ParallelBintScoreInitialPoseAndTrajectory("TrajBintsCube", title="Trajectory Binding Interactions")
    clusCube = ParallelClusterOETrajCube("ClusterOETrajCube", title="Clustering")
    clusPop = ParallelClusterPopAnalysis('ClusterPopAnalysis', title="Clustering Analysis")
    clusOEMols = ParallelMakeClusterTrajOEMols('MakeClusterTrajOEMols', title="Per-Cluster Analysis")
    prepDataset = ParallelTrajAnalysisReportDataset('TrajAnalysisReportDataset', title="Analysis Report")
    report_gen = ParallelMDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport", title="Relevant Output Extraction")
    analysis_group = ParallelCubeGroup(cubes=[catLigTraj, catLigMMPBSA, trajBints, clusCube, clusPop,
                                              clusOEMols, prepDataset, report_gen])
    input_floe.add_group(analysis_group)

    report = MDFloeReportCube("report", title="Floe Report")

    input_floe.add_cubes(trajCube, IntECube, PBSACube, confGather,
                         catLigTraj, catLigMMPBSA, trajBints, clusCube, clusPop, clusOEMols,
                         prepDataset, report_gen, report)
    # Success Connections
    input_cube.success.connect(trajCube.intake)
    trajCube.success.connect(IntECube.intake)
    IntECube.success.connect(PBSACube.intake)
    PBSACube.success.connect(confGather.intake)
    confGather.success.connect(catLigTraj.intake)
    catLigTraj.success.connect(catLigMMPBSA.intake)
    catLigMMPBSA.success.connect(trajBints.intake)
    trajBints.success.connect(clusCube.intake)
    clusCube.success.connect(clusPop.intake)
    clusPop.success.connect(clusOEMols.intake)
    clusOEMols.success.connect(prepDataset.intake)
    prepDataset.success.connect(report_gen.intake)
    report_gen.success.connect(report.intake)
    report.du_success.connect(du_output_cube.intake)

    # Fail Connections
    trajCube.failure.connect(fail_cube.fail_in)
    IntECube.failure.connect(fail_cube.fail_in)
    PBSACube.failure.connect(fail_cube.fail_in)
    confGather.failure.connect(fail_cube.fail_in)
    catLigTraj.failure.connect(fail_cube.fail_in)
    catLigMMPBSA.failure.connect(fail_cube.fail_in)
    trajBints.failure.connect(fail_cube.fail_in)
    clusCube.failure.connect(fail_cube.fail_in)
    clusPop.failure.connect(fail_cube.fail_in)
    clusOEMols.failure.connect(fail_cube.fail_in)
    prepDataset.failure.connect(fail_cube.fail_in)
    report_gen.failure.connect(fail_cube.fail_in)

    return report

