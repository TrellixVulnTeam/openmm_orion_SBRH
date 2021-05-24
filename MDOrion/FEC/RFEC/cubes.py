# (C) 2020 OpenEye Scientific Software Inc. All rights reserved.
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


from orionplatform.mixins import RecordPortsMixin

from orionplatform.ports import (RecordInputPort,
                                 RecordOutputPort)

from floe.api import (parameters,
                      ComputeCube,
                      ParallelMixin)

from floereport import FloeReport, LocalFloeReport

import MDOrion.TrjAnalysis.Affinity_FloeReport_utils as flrpt

from MDOrion.FEC.RFEC.utils import CombineAndOffsetPredExptDGs

from MDOrion.Standards.standards import Fields

import traceback

from orionplatform.parameters import FileInputParameter

from orionclient.utils import TemporaryPath

from MDOrion.FEC.RFEC import utils

from datarecord import OERecord

from oemdtoolbox.FEC.RBFEC.chimera import Chimera

from MDOrion.FEC.RFEC.utils import (gmx_chimera_topology_injection,
                                    gmx_chimera_coordinate_injection,
                                    upload_gmx_files,
                                    download_gmx_file)

from MDOrion.FEC.RFEC.utils import parmed_find_ligand


from datarecord import (Types,
                        OEField,
                        Meta,
                        OEFieldMeta)


from MDOrion.Standards.mdrecord import MDDataRecord

from MDOrion.MDEngines.utils import MDState

import numpy as np

import os

from os import environ

from simtk import unit

import io

from MDOrion.FEC.RFEC.gmx_run import check_gmx_grompp

import tempfile

import tarfile

from orionclient.session import in_orion


class BoundUnboundSwitchCube(RecordPortsMixin, ComputeCube):
    title = "Bound and UnBound Switching Cube"
    
    classification = [["Simulation Flask Preparation"]]
    tags = ['Simulation', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube emits complexes on the bound port and non-complexes on
    the standard out port. The Flask ids are re-set
    """

    uuid = "80e656e8-33c5-4560-99b5-95e4f69c1701"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    bound_port = RecordOutputPort("bound_port", initializer=False)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count = 0

    def process(self, record, port):

        try:
            if not record.has_value(Fields.md_components):
                raise ValueError("MD Components Field is missing")

            md_components = record.get_value(Fields.md_components)

            if not record.has_value(Fields.flaskid):
                record.set_value(Fields.flaskid, self.count)
                self.count += 1

            if md_components.has_protein:
                if not record.has_value(Fields.FEC.RBFEC.thd_leg_type):
                    record.set_value(Fields.FEC.RBFEC.thd_leg_type, "Bound_OPLMD")
                self.bound_port.emit(record)
            else:
                if not record.has_value(Fields.FEC.RBFEC.thd_leg_type):
                    record.set_value(Fields.FEC.RBFEC.thd_leg_type, "UnBound_OPLMD")
                self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class RBFECEdgeGathering(RecordPortsMixin, ComputeCube):
    title = "RBFEC Edge Gathering"

    classification = [["Relative Free Energy"]]
    tags = ['Ligand', 'Edge Mapping']
    description = """
    This cube gathers the MD equilibrium simulation data related to the Bound and Unbound simulations
    and creates edges used to perform relative binding affinity calculations. Each edge is created
    providing an edge mapping text file where each row contains a string with the ligand 
    names involved in the relative binding free energy calculation following the syntax:
    ....
    lig_i_name >> lig_j_name
    ....
    The cube produces as output records with the Bound and Unbound edge data

    Input:
    -------
    Data record Stream - Streamed-in of the Bound simulation data
    Data record Stream - Streamed-in of the Unbound simulation data
    File name - File with the relative binding free energy mapping
 
    Output:
    -------
    Data Record Stream - Streamed-out of records where Bound and Unbound 
    data has been paired to run relative binding free energy calculations.
    """

    uuid = "970052ea-25dc-4aa6-87ee-840506022a8b"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    map_file = FileInputParameter("map_file", title="RBFEC Mapping file",
                                  description="RBFEC mapping file", required=True,
                                  default=None)

    bound_port = RecordInputPort("bound_port", initializer=False)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

        self.Bound_edges = dict()
        self.Bound_records = dict()

        self.Unbound_edges = dict()
        self.Unbound_records = dict()

        self.Bound_Unbound_edges = dict()

        file = list(self.args.map_file)
        try:
            for file_obj in file:
                with TemporaryPath() as path:
                    file_obj.copy_to(path)
                    with open(path, "r") as f:
                        edge_list = f.readlines()
        except Exception as e:
            raise IOError("Problems detected reading the edge mapping file: {}".format(str(e)))

        if not edge_list:
            raise IOError("Edge file is empty {}")

        for edge in edge_list:

            # Empty line
            if not edge:
                continue

            # Comment
            if edge.startswith(";") or edge.startswith("#"):
                continue

            # New line
            if edge == '\n':
                continue

            edge_str = utils.edge_map_grammar(edge)

            if len(edge_str) > 3:
                raise ValueError("Syntax Error Edge File: {}".format(edge_str))

            edge_A_State = edge_str[0]
            edge_B_State = edge_str[2]

            if edge_A_State == edge_B_State:
                self.opt['Logger'].warn("Edge with the same starting and final state detected. "
                                        "The edge will be skipped")
                continue

            self.count = 0
            self.Bound_edges[(edge_A_State, edge_B_State)] = [False, False]
            self.Unbound_edges[(edge_A_State, edge_B_State)] = [False, False]
            self.Bound_Unbound_edges[(edge_A_State, edge_B_State)] = [self.Bound_edges, self.Unbound_edges]

        if not self.Bound_edges:
            raise ValueError("Not available edges have been found. Please check your edge file")

    def process(self, record, port):
        try:

            if port == 'intake':
                edges = self.Unbound_edges
                lig_records = self.Unbound_records

            else:
                edges = self.Bound_edges
                lig_records = self.Bound_records

            if not record.has_value(Fields.ligand_name):
                raise ValueError("The record is missing the ligand name field")

            lig_name = record.get_value(Fields.ligand_name)

            # Populate the ligand Bound/Unbound record dictionary
            lig_records[lig_name] = record

            # Mark True all the occurrences of lig_name in the A and B state
            for edge, occurrence in edges.items():

                lig_A_name = edge[0]
                lig_B_name = edge[1]

                if lig_name == lig_A_name:
                    occurrence[0] = True

                elif lig_name == lig_B_name:
                    occurrence[1] = True

            full_edges_to_del = [(edge[0], edge[1]) for edge, v in self.Bound_Unbound_edges.items() if v[0][edge][0]
                                 and v[0][edge][1] and v[1][edge][0] and v[1][edge][1]]

            # full_edges_to_del = []
            # for ed, v in self.Bound_Unbound_edges.items():
            #     if v[0][ed][0] and v[0][ed][1] and v[1][ed][0] and v[1][ed][1]:
            #         full_edges_to_del.append((ed[0], ed[1]))

            for edge_to_del in full_edges_to_del:
                if edge_to_del in edges:

                    del self.Bound_edges[edge_to_del]
                    del self.Unbound_edges[edge_to_del]
                    del self.Bound_Unbound_edges[edge_to_del]

                    lig_A_name = edge_to_del[0]
                    lig_B_name = edge_to_del[1]

                    bound_rec_A = self.Bound_records[lig_A_name]
                    bound_rec_B = self.Bound_records[lig_B_name]

                    unbound_rec_A = self.Unbound_records[lig_A_name]
                    unbound_rec_B = self.Unbound_records[lig_B_name]

                    new_record = OERecord()
                    new_record.set_value(Fields.FEC.RBFEC.NESC.state_A, [bound_rec_A, unbound_rec_A])
                    new_record.set_value(Fields.FEC.RBFEC.NESC.state_B, [bound_rec_B, unbound_rec_B])
                    new_record.set_value(Fields.FEC.RBFEC.edge_name, lig_A_name + '_to_' + lig_B_name)
                    new_record.set_value(Fields.FEC.RBFEC.edgeid, self.count)
                    self.success.emit(new_record)
                    self.count += 1

            # Clean Bound and Unbound Records
            bound_rec_to_del = list()
            for ln in self.Bound_records.keys():
                to_del = True
                for ed in self.Bound_edges:
                    if ln == ed[0] or ln == ed[1]:
                        to_del = False
                if to_del:
                    bound_rec_to_del.append(ln)

            for ln in bound_rec_to_del:
                del self.Bound_records[ln]

            unbound_rec_to_del = list()
            for ln in self.Unbound_records.keys():
                to_del = True
                for ed in self.Unbound_edges:
                    if ln == ed[0] or ln == ed[1]:
                        to_del = False
                if to_del:
                    unbound_rec_to_del.append(ln)

            for ln in unbound_rec_to_del:
                del self.Unbound_records[ln]

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):

        if self.Bound_edges:
            self.opt['Logger'].warn("The Bound edge list is not empty :{}".format(self.Bound_edges))

        if self.Unbound_edges:
            self.opt['Logger'].warn("The UnBound edge list is not empty :{}".format(self.Unbound_edges))

        return


class NESGMXChimera(RecordPortsMixin, ComputeCube):
    title = "GMX Chimera"

    classification = [["Relative Free Energy"]]
    tags = ['Ligand', 'Edge Mapping']
    description = """
    Relative Binding Free Energy calculations via Alchemical Transformations are usually
    done by using a chimeric molecule where the common scaffold between two ligands
    involved in an edge is modified to contain the atom excesses of the two edge ligands. This
    cube generates the chimera molecule and inject it on select frames from the Bound and 
    Unbound simulations for the Gromacs forward and reverse Non Equilibrium simulations.  
    
    Input:
    -------
    Data record Stream - Streamed-in of the Gathered Bound and Unbound simulation data
   
    Output:
    -------
    Data Record Stream - Streamed-out of records containing Gromacs forward and reverse
    Bound and Unbound data.
    """

    uuid = "0ccfad69-bada-48fa-a890-348d7784d7ea"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 32000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    # Ligand Residue Name
    lig_res_name = parameters.StringParameter('lig_res_name',
                                              default='LIG',
                                              help_text='The new ligand residue name')

    trajectory_frames = parameters.IntegerParameter(
        'trajectory_frames',
        default=80,
        help_text="The total number of trajectory frames to run NES")

    bound_port = RecordOutputPort("bound_port", initializer=False)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:

            edge_id = record.get_value(Fields.FEC.RBFEC.edgeid)
            edge_name = record.get_value(Fields.FEC.RBFEC.edge_name)

            rec_list_state_A = record.get_value(Fields.FEC.RBFEC.NESC.state_A)
            rec_list_state_B = record.get_value(Fields.FEC.RBFEC.NESC.state_B)

            self.opt['Logger'].info("[{}] Processing Edge: {}".format(self.title, edge_name))

            md_record_state_A_Bound = MDDataRecord(rec_list_state_A[0])
            md_record_state_B_Bound = MDDataRecord(rec_list_state_B[0])

            md_record_state_A_Unbound = MDDataRecord(rec_list_state_A[1])
            md_record_state_B_Unbound = MDDataRecord(rec_list_state_B[1])

            lig_A = md_record_state_A_Unbound.get_md_components.get_ligand
            lig_B = md_record_state_B_Unbound.get_md_components.get_ligand

            pmd_flask_state_A_Unbound = md_record_state_A_Unbound.get_parmed(sync_stage_name="Flask Parametrization")
            pmd_flask_state_B_Unbound = md_record_state_B_Unbound.get_parmed(sync_stage_name="Flask Parametrization")

            pmd_flask_state_A_Bound = md_record_state_A_Bound.get_parmed(sync_stage_name="Flask Parametrization")
            pmd_flask_state_B_Bound = md_record_state_B_Bound.get_parmed(sync_stage_name="Flask Parametrization")

            pmd_lig_A, idxA = parmed_find_ligand(pmd_flask_state_A_Unbound)
            pmd_lig_B, idxB = parmed_find_ligand(pmd_flask_state_B_Unbound)

            if pmd_lig_A is None:
                raise ValueError("It was not possible to extract the ligand parmed from the state A")

            if pmd_lig_B is None:
                raise ValueError("It was not possible to extract the ligand parmed from the state B")

            chimera = Chimera(lig_A, lig_B, pmd_lig_A, pmd_lig_B)

            pmd_chimera_A_to_B_initial, pmd_chimera_A_to_B_final, graph_A_to_B_dic = chimera.pmd_chimera(morph="A_to_B")
            pmd_chimera_B_to_A_initial, pmd_chimera_B_to_A_final, graph_B_to_A_dic = chimera.pmd_chimera(morph="B_to_A")

            gmx_top_A_to_B_Unbound = gmx_chimera_topology_injection(pmd_flask_state_A_Unbound,
                                                                    pmd_chimera_A_to_B_initial,
                                                                    pmd_chimera_A_to_B_final)

            gmx_top_B_to_A_Unbound = gmx_chimera_topology_injection(pmd_flask_state_B_Unbound,
                                                                    pmd_chimera_B_to_A_initial,
                                                                    pmd_chimera_B_to_A_final)

            gmx_top_A_to_B_Bound = gmx_chimera_topology_injection(pmd_flask_state_A_Bound,
                                                                  pmd_chimera_A_to_B_initial,
                                                                  pmd_chimera_A_to_B_final)

            gmx_top_B_to_A_Bound = gmx_chimera_topology_injection(pmd_flask_state_B_Bound,
                                                                  pmd_chimera_B_to_A_initial,
                                                                  pmd_chimera_B_to_A_final)

            gmx_gro_A_to_B_Unbound = gmx_chimera_coordinate_injection(pmd_chimera_A_to_B_initial,
                                                                      md_record_state_A_Unbound,
                                                                      self.opt['trajectory_frames'],
                                                                      lig_B,
                                                                      chimera,
                                                                      graph_A_to_B_dic)

            gmx_gro_B_to_A_Unbound = gmx_chimera_coordinate_injection(pmd_chimera_B_to_A_initial,
                                                                      md_record_state_B_Unbound,
                                                                      self.opt['trajectory_frames'],
                                                                      lig_A,
                                                                      chimera,
                                                                      graph_B_to_A_dic)

            gmx_gro_A_to_B_Bound = gmx_chimera_coordinate_injection(pmd_chimera_A_to_B_initial,
                                                                    md_record_state_A_Bound,
                                                                    self.opt['trajectory_frames'],
                                                                    lig_B,
                                                                    chimera,
                                                                    graph_A_to_B_dic)

            gmx_gro_B_to_A_Bound = gmx_chimera_coordinate_injection(pmd_chimera_B_to_A_initial,
                                                                    md_record_state_B_Bound,
                                                                    self.opt['trajectory_frames'],
                                                                    lig_A,
                                                                    chimera,
                                                                    graph_B_to_A_dic)

            self.opt['Logger'].info("GMX Chimera Topology Created")

            frame_count_field = OEField("frame_count", Types.Int)

            for count, gmx_gro in enumerate(gmx_gro_A_to_B_Unbound):

                if count == 0:
                    check_gmx_grompp(gmx_gro, gmx_top_A_to_B_Unbound, sim_type=edge_name + ' Forward Unbound')

                if in_orion():
                    gmx_tar_A_to_B_Unbound = tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="gmx_A_to_B_un_", suffix='_'+str(count)+".tar")
                else:
                    gmx_tar_A_to_B_Unbound = tempfile.NamedTemporaryFile(mode='w', dir="./", delete=False, prefix="gmx_A_to_B_un_", suffix='_' + str(count) + ".tar")

                with tarfile.open(gmx_tar_A_to_B_Unbound.name, mode='w:gz') as archive:

                    with tempfile.TemporaryDirectory() as outdir:
                        gmx_top_fn = os.path.join(outdir, "gmx_top.top")
                        gmx_gro_fn = os.path.join(outdir, "gmx_gro.gro")

                        with open(gmx_top_fn, 'w') as f:
                            f.write(gmx_top_A_to_B_Unbound)
                        with open(gmx_gro_fn, 'w') as f:
                            f.write(gmx_gro)

                        archive.add(gmx_top_fn, arcname=os.path.basename(gmx_top_fn))
                        archive.add(gmx_gro_fn, arcname=os.path.basename(gmx_gro_fn))

                if in_orion():
                    if not upload_gmx_files(gmx_tar_A_to_B_Unbound.name, md_record_state_A_Unbound, shard_name=os.path.basename(gmx_tar_A_to_B_Unbound.name)):
                        raise ValueError("It was not possible to upload the gromacs file to Orion")
                else:
                    if not upload_gmx_files(os.path.basename(gmx_tar_A_to_B_Unbound.name), md_record_state_A_Unbound, shard_name=os.path.basename(gmx_tar_A_to_B_Unbound.name)):
                        raise ValueError("It was not possible to set the gromacs file")

                gmx_tar_A_to_B_Unbound.close()

                md_record_state_A_Unbound.set_value(Fields.FEC.RBFEC.edgeid, edge_id)
                md_record_state_A_Unbound.set_value(Fields.FEC.RBFEC.edge_name, edge_name)
                md_record_state_A_Unbound.set_value(Fields.FEC.RBFEC.NESC.direction, "Forward_OPLMD")
                md_record_state_A_Unbound.set_value(frame_count_field, count)

                self.success.emit(md_record_state_A_Unbound.get_record)

            self.opt['Logger'].info("GMX Chimera Unbound {} Forward Coordinate Created".format(edge_name))

            for count, gmx_gro in enumerate(gmx_gro_A_to_B_Bound):

                if count == 0:
                    check_gmx_grompp(gmx_gro, gmx_top_A_to_B_Bound, sim_type=edge_name + ' Forward Bound')

                if in_orion():
                    gmx_tar_A_to_B_Bound = tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="gmx_A_to_B_bn_", suffix='_' + str(count) + ".tar")
                else:
                    gmx_tar_A_to_B_Bound = tempfile.NamedTemporaryFile(mode='w', dir="./", delete=False, prefix="gmx_A_to_B_bn_", suffix='_' + str(count) + ".tar")

                with tarfile.open(gmx_tar_A_to_B_Bound.name, mode='w:gz') as archive:

                    with tempfile.TemporaryDirectory() as outdir:
                        gmx_top_fn = os.path.join(outdir, "gmx_top.top")
                        gmx_gro_fn = os.path.join(outdir, "gmx_gro.gro")

                        with open(gmx_top_fn, 'w') as f:
                            f.write(gmx_top_A_to_B_Bound)
                        with open(gmx_gro_fn, 'w') as f:
                            f.write(gmx_gro)

                        archive.add(gmx_top_fn, arcname=os.path.basename(gmx_top_fn))
                        archive.add(gmx_gro_fn, arcname=os.path.basename(gmx_gro_fn))

                if in_orion():
                    if not upload_gmx_files(gmx_tar_A_to_B_Bound.name, md_record_state_A_Bound, shard_name=os.path.basename(gmx_tar_A_to_B_Bound.name)):
                        raise ValueError("It was not possible to upload the gromacs file to Orion")

                else:
                    if not upload_gmx_files(os.path.basename(gmx_tar_A_to_B_Bound.name), md_record_state_A_Bound, shard_name=os.path.basename(gmx_tar_A_to_B_Bound.name)):
                        raise ValueError("It was not possible to set the gromacs file to Orion")

                gmx_tar_A_to_B_Bound.close()

                md_record_state_A_Bound.set_value(Fields.FEC.RBFEC.edgeid, edge_id)
                md_record_state_A_Bound.set_value(Fields.FEC.RBFEC.edge_name, edge_name)
                md_record_state_A_Bound.set_value(Fields.FEC.RBFEC.NESC.direction, "Forward_OPLMD")
                md_record_state_A_Bound.set_value(frame_count_field, count)
                self.bound_port.emit(md_record_state_A_Bound.get_record)

            self.opt['Logger'].info("GMX Chimera Bound {} Forward Coordinate Created".format(edge_name))

            b_to_a_name = edge_name.split("_to_")[1] + "_to_" + edge_name.split("_to_")[0]

            for count, gmx_gro in enumerate(gmx_gro_B_to_A_Unbound):

                if count == 0:
                    check_gmx_grompp(gmx_gro, gmx_top_B_to_A_Unbound, sim_type=edge_name + ' Reverse Unbound')

                if in_orion():
                    gmx_tar_B_to_A_Unbound = tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="gmx_B_to_A_un_", suffix='_' + str(count) + ".tar")

                else:
                    gmx_tar_B_to_A_Unbound = tempfile.NamedTemporaryFile(mode='w', dir="./", delete=False, prefix="gmx_B_to_A_un_", suffix='_' + str(count) + ".tar")

                with tarfile.open(gmx_tar_B_to_A_Unbound.name, mode='w:gz') as archive:

                    with tempfile.TemporaryDirectory() as outdir:
                        gmx_top_fn = os.path.join(outdir, "gmx_top.top")
                        gmx_gro_fn = os.path.join(outdir, "gmx_gro.gro")

                        with open(gmx_top_fn, 'w') as f:
                            f.write(gmx_top_B_to_A_Unbound)
                        with open(gmx_gro_fn, 'w') as f:
                            f.write(gmx_gro)

                        archive.add(gmx_top_fn, arcname=os.path.basename(gmx_top_fn))
                        archive.add(gmx_gro_fn, arcname=os.path.basename(gmx_gro_fn))

                if in_orion():
                    if not upload_gmx_files(gmx_tar_B_to_A_Unbound.name, md_record_state_B_Unbound, shard_name=os.path.basename(gmx_tar_B_to_A_Unbound.name)):
                        raise ValueError("It was not possible to upload the gromacs file to Orion")

                else:
                    if not upload_gmx_files(os.path.basename(gmx_tar_B_to_A_Unbound.name), md_record_state_B_Unbound, shard_name=os.path.basename(gmx_tar_B_to_A_Unbound.name)):
                        raise ValueError("It was not possible to set the gromacs file to Orion")

                gmx_tar_B_to_A_Unbound.close()

                md_record_state_B_Unbound.set_value(Fields.FEC.RBFEC.edgeid, edge_id)
                md_record_state_B_Unbound.set_value(Fields.FEC.RBFEC.edge_name, b_to_a_name)
                md_record_state_B_Unbound.set_value(Fields.FEC.RBFEC.NESC.direction, "Reverse_OPLMD")
                md_record_state_B_Unbound.set_value(frame_count_field, count)
                self.success.emit(md_record_state_B_Unbound.get_record)

            self.opt['Logger'].info("GMX Chimera Unbound {} Reverse Coordinate Created".format(edge_name))

            for count, gmx_gro in enumerate(gmx_gro_B_to_A_Bound):

                if count == 0:
                    check_gmx_grompp(gmx_gro, gmx_top_B_to_A_Bound, sim_type=edge_name + ' Reverse Bound')

                if in_orion():
                    gmx_tar_B_to_A_Bound = tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="gmx_B_to_A_bn_", suffix='_' + str(count) + ".tar")
                else:
                    gmx_tar_B_to_A_Bound = tempfile.NamedTemporaryFile(mode='w', dir="./", delete=False, prefix="gmx_B_to_A_bn_", suffix='_' + str(count) + ".tar")
                with tarfile.open(gmx_tar_B_to_A_Bound.name, mode='w:gz') as archive:

                    with tempfile.TemporaryDirectory() as outdir:
                        gmx_top_fn = os.path.join(outdir, "gmx_top.top")
                        gmx_gro_fn = os.path.join(outdir, "gmx_gro.gro")

                        with open(gmx_top_fn, 'w') as f:
                            f.write(gmx_top_B_to_A_Bound)
                        with open(gmx_gro_fn, 'w') as f:
                            f.write(gmx_gro)

                        archive.add(gmx_top_fn, arcname=os.path.basename(gmx_top_fn))
                        archive.add(gmx_gro_fn, arcname=os.path.basename(gmx_gro_fn))

                if in_orion():
                    if not upload_gmx_files(gmx_tar_B_to_A_Bound.name, md_record_state_B_Bound, shard_name=os.path.basename(gmx_tar_B_to_A_Bound.name)):
                        raise ValueError("It was not possible to upload the gromacs file to Orion")
                else:
                    if not upload_gmx_files(os.path.basename(gmx_tar_B_to_A_Bound.name), md_record_state_B_Bound, shard_name=os.path.basename(gmx_tar_B_to_A_Bound.name)):
                        raise ValueError("It was not possible to set the gromacs file to Orion")

                gmx_tar_B_to_A_Bound.close()

                md_record_state_B_Bound.set_value(Fields.FEC.RBFEC.edgeid, edge_id)
                md_record_state_B_Bound.set_value(Fields.FEC.RBFEC.edge_name, b_to_a_name)
                md_record_state_B_Bound.set_value(Fields.FEC.RBFEC.NESC.direction, "Reverse_OPLMD")
                md_record_state_B_Bound.set_value(frame_count_field, count)
                self.bound_port.emit(md_record_state_B_Bound.get_record)

            self.opt['Logger'].info("GMX Chimera Bound {} Reverse Coordinate Created".format(edge_name))
            self.opt['Logger'].info("GMX Chimera {} Processed".format(edge_name))

        except Exception as e:
            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)


class NESGMX(RecordPortsMixin, ComputeCube):
    title = "NES GMX"
    
    classification = [["Free Energy"]]
    tags = ["Ligand", "Protein", "Free Energy", "Non Equilibrium"]
    description = """
    This cubes runs the forward and reverse Bound and Unbound starting equilibrium
    frames by using Gromacs for each edge. The physical work associated with the 
    edge mutations are calculated by using the Thermodynamic Integration method 
    and saved on the record ready to be analyzed.
    
    Input:
    -------
    Data record Stream - Streamed-in of records containing Gromacs equilibrium frames 
    data for the forward and reverse Bound and Unbound calculations 
   
    Output:
    -------
    Data Record Stream - Streamed-out of records containing work data used to generate relative
    binding affinity data.
    """

    uuid = "3641fe19-780f-4998-90c5-2ec4102121ba"

    # Override defaults for some parameters

    # parameter_overrides = {
    #     "gpu_count": {"default": 1},
    #     "instance_type": {"default": "g3.4xlarge"},  # Gpu Family selection
    #     "memory_mb": {"default": 14000},
    #     "spot_policy": {"default": "Required"},
    #     "prefetch_count": {"default": 1},  # 1 molecule at a time
    #     "item_count": {"default": 1},  # 1 molecule at a time
    #     "max_failures": {"default": 2}  # it is going to retry just one more time
    # }

    parameter_overrides = {
        "cpu_count": {"default": 16},
        "memory_mb": {"default": float(8 * 1.8 * 1024)},
        "gpu_count": {"default": 1},
        "disk_space": {"default": float(6.0 * 1024)},
        "spot_policy": {"default": "Required"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 2}  # just one retry
    }

    temperature = parameters.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameters.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    time = parameters.DecimalParameter(
        'time',
        default=0.05,
        help_text="NPT simulation time in nanoseconds")

    constraints = parameters.StringParameter(
        'constraints',
        default='All-Bonds',
        choices=['None', 'Bonds2H', 'Angles2H', 'All-Bonds'],
        help_text="""None, Bonds2H, Angles2H, or All-Bonds
            Which type of constraints to add to the system.
            None means no bonds are constrained.
            Bonds2H means bonds with hydrogen are constrained, etc.""")

    restraints = parameters.BooleanParameter(
        'restraints',
        default=True,
        help_text=""""If True restraints are applied to 
        the chimeric molecule during the gromacs equilibration""")

    verbose = parameters.BooleanParameter(
        'verbose',
        default=False,
        help_text='Increase log file verbosity')

    suffix = parameters.StringParameter(
        'suffix',
        default='nes',
        help_text='Filename suffix for output simulation files')

    gmx_process_timeout = parameters.DecimalParameter(
        'gmx_process_timeout',
        default=3600.0,
        help_text="Gromacs process timeout in seconds")

    gmx_openmp_threads = parameters.IntegerParameter(
        'gmx_openmp_threads',
        default=16,
        help_text='Number of Gromacs OpenMP threads')

    gmx_mpi_threads = parameters.IntegerParameter(
        'gmx_mpi_threads',
        default=1,
        help_text='Number of Gromacs MPI threads')

    def begin(self):
            self.opt = vars(self.args)
            self.opt['Logger'] = self.log
            self.edge_dic = dict()

    def process(self, record, port):

        try:

            opt = dict(self.opt)
            opt['CubeTitle'] = self.title

            if not record.has_field(Fields.title):
                raise ValueError("Missing title field")

            flask_title = record.get_value(Fields.title)

            frame_count = record.get_value(OEField("frame_count", Types.Int))

            mdrecord = MDDataRecord(record)

            extra_data_fn = download_gmx_file(mdrecord)

            with tarfile.open(extra_data_fn) as tar:
                tar.extractall(path=mdrecord.cwd)

            gmx_gro_fn = os.path.join(mdrecord.cwd, "gmx_gro.gro")
            with open(gmx_gro_fn, 'r') as f:
                gmx_gro_str = f.read()

            gmx_top_fn = os.path.join(mdrecord.cwd, "gmx_top.top")
            with open(gmx_top_fn, 'r') as f:
                gmx_top_str = f.read()

            md_components = mdrecord.get_md_components

            opt['frame_count'] = frame_count
            opt['out_directory'] = mdrecord.cwd
            opt['out_prefix'] = os.path.basename(mdrecord.cwd)+'_'+flask_title+'_'+str(frame_count)
            opt['trj_fn'] = opt['out_prefix'] + '_' + opt['suffix'] + '_' + 'traj.tar.gz'
            # TODO This is not used for now. NES Trajectories are not uploaded
            trj_fn = opt['trj_fn']

            mdstate = mdrecord.get_stage_state(stg_name='last')

            box = mdstate.get_box_vectors()

            if md_components.get_box_vectors is not None:

                box_v = box.value_in_unit(unit.angstrom)
                box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])

                min_box = np.min(box_v)
                opt['min_box'] = min_box

            # Run Gromacs
            utils.gmx_nes_run(gmx_gro_str, gmx_top_str, opt)

            str_logger = '\n' + '-' * 32 + ' SIMULATION FEC NE' + '-' * 32

            with(io.open(os.path.join(opt['out_directory'], opt['log_fn']), 'r', encoding='utf8', errors='ignore')) as flog:
                str_logger += '\n' + flog.read()

            data_fn = opt['out_prefix'] + '.tar.gz'

            # The Parmed structure, the flask, the gromacs positions and the work
            # are updated inside the Gromacs NES MD run
            mdstate = MDState(opt['pmd'])

            # TODO DISABLED FOR NOW
            # if not mdrecord.add_new_stage(self.title,
            #                               MDStageTypes.FEC,
            #                               opt['flask'],
            #                               mdstate,
            #                               data_fn,
            #                               append=True,
            #                               log=str_logger,
            #                               # trajectory_fn=trj_fn,
            #                               # trajectory_engine=MDEngines.Gromacs,
            #                               # trajectory_orion_ui=flask_title + '_' + str(frame_count)
            #                               ):
            #
            #     raise ValueError("Problems adding in the new FEC Stage")
            #
            # # Update Gromacs coordinates on the record
            # record.set_value(Fields.FEC.RBFEC.NESC.gmx_gro, opt['gro_str'])

            # Set the calculated work
            record.set_value(Fields.FEC.RBFEC.NESC.work, opt['gmx_work'])

            self.success.emit(mdrecord.get_record)

            # del mdrecord

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class NESAnalysis(RecordPortsMixin, ComputeCube):
    title = "NES Analysis"
    classification = [["FEC Analysis"]]
    tags = ['Complex', 'Protein', 'Ligand', 'Solvation']
    description = """
    The Analysis cube collects the data generated from the NES runs and calculates
    the relative binding affinities for the edges involved in the Relative
    binding affinity calculations. The binding affinities are estimated by using 
    the Bennet Acceptance Ratio method.
    
    Input:
    -------
    Data Record Stream - Streamed-in of records containing work data produced by NES runs
     
    Output:
    -------
    Data record Stream - Streamed-out of records containing relative binding affinity data 
    for each edge.
   
    """

    uuid = "50ccc16d-67ae-4b4f-9a98-2e6b8ecb1868"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 16000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameters.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    units = parameters.StringParameter(
        'units',
        choices=['kcal/mol', 'kJ/mol'],
        default='kcal/mol',
        help_text='Units to use to display the plots'
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.edgeid_works = dict()
        self.edgeid_ligands = dict()
        self.edgeid_lig_names = dict()
        self.collections = dict()

    def process(self, record, port):

        try:
            lig_name = record.get_value(Fields.ligand_name)
            md_components = record.get_value(Fields.md_components)
            ligand = md_components.get_ligand
            ligand.SetTitle(lig_name)
            leg_type = record.get_value(Fields.FEC.RBFEC.thd_leg_type)
            edgeid = record.get_value(Fields.FEC.RBFEC.edgeid)

            if not len(self.collections):
                self.collections = record.get_value(Fields.collections)

            if edgeid not in self.edgeid_works:
                work_forward_bound = dict()
                work_reverse_bound = dict()
                work_forward_unbound = dict()
                work_reverse_unbound = dict()

                self.edgeid_works[edgeid] = [work_forward_bound,
                                             work_reverse_bound,
                                             work_forward_unbound,
                                             work_reverse_unbound]

                self.edgeid_ligands[edgeid] = [None, None]
                self.edgeid_lig_names[edgeid] = [None, None]

            frame_count = record.get_value(Fields.FEC.RBFEC.NESC.frame_count)
            work = record.get_value(Fields.FEC.RBFEC.NESC.work)

            # TODO change the OEField in direction
            direction = record.get_value(Fields.FEC.RBFEC.NESC.direction)

            # Forward
            if direction == 'Forward_OPLMD':

                # Both ligands in the edge have been already collected
                if self.edgeid_ligands[edgeid][0] and self.edgeid_ligands[edgeid][1]:
                    pass
                elif self.edgeid_ligands[edgeid][0] is None:
                    self.edgeid_ligands[edgeid][0] = ligand
                    self.edgeid_lig_names[edgeid][0] = lig_name
                if leg_type == "Bound_OPLMD":
                    self.edgeid_works[edgeid][0][frame_count] = work
                else:
                    self.edgeid_works[edgeid][2][frame_count] = work

            # Reverse
            else:

                if self.edgeid_ligands[edgeid][0] and self.edgeid_ligands[edgeid][1]:
                    pass
                elif self.edgeid_ligands[edgeid][1] is None:
                    self.edgeid_ligands[edgeid][1] = ligand
                    self.edgeid_lig_names[edgeid][1] = lig_name

                if leg_type == "Bound_OPLMD":
                    self.edgeid_works[edgeid][1][frame_count] = work
                else:
                    self.edgeid_works[edgeid][3][frame_count] = work

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):
        count = 0

        try:
            for edgeid, work_list in self.edgeid_works.items():
                try:

                    if self.edgeid_ligands[edgeid][0] is None and self.edgeid_ligands[edgeid][1] is None:
                        self.opt['Logger'].warn(
                            "The edge id {} is missing ligand data in the A and B States. "
                            "The edge will be skipped: {}".format(edgeid, self.edgeid_lig_names[edgeid]),
                            flush=True)
                        continue
                    elif self.edgeid_ligands[edgeid][0] is None and self.edgeid_ligands[edgeid][1]:
                        self.opt['Logger'].warn(
                            "The edge id {} - (None, {}) is missing ligand data in the A State. "
                            "The edge will be skipped: {}".format(edgeid, self.edgeid_lig_names[edgeid][1],
                                                                  self.edgeid_lig_names[edgeid]),
                            flush=True)
                        continue
                    elif self.edgeid_ligands[edgeid][0] and self.edgeid_ligands[edgeid][1] is None:
                        self.opt['Logger'].warn(
                            "The edge id {} - ({}, None) is missing ligand data in the B State. "
                            "The edge will be skipped: {}".format(edgeid, self.edgeid_lig_names[edgeid][0],
                                                                  self.edgeid_lig_names[edgeid]),
                            flush=True)
                        continue

                    self.opt['Logger'].info("...Processing Edge {} - {}".
                                            format(edgeid, self.edgeid_lig_names[edgeid]), flush=True)

                    # Make edge depiction for the floe report tile
                    ligandA = self.edgeid_ligands[edgeid][0]
                    ligandB = self.edgeid_ligands[edgeid][1]

                    # The work dictionaries are ordered by frame counts and the reverse
                    # works are sign inverted. Work values in kJ/mol

                    if None in work_list:
                        self.opt['Logger'].warn('The edge id {} is missing work data and it will be skipped'.
                                                format(edgeid),  flush=True)
                        continue

                    # Ordered Dictionary frame:work
                    forward_bound = {k: work_list[0][k] for k in sorted(work_list[0])}
                    reverse_bound = {k: -work_list[1][k] for k in sorted(work_list[1])}
                    forward_unbound = {k: work_list[2][k] for k in sorted(work_list[2])}
                    reverse_unbound = {k: -work_list[3][k] for k in sorted(work_list[3])}

                    # NES Data Analysis. Mute output
                    # with open(os.devnull, 'w') as devnull:
                    #     with contextlib.redirect_stdout(devnull):
                    #         results = utils.nes_data_analysis(forward_bound,
                    #                                           reverse_bound,
                    #                                           forward_unbound,
                    #                                           reverse_unbound)

                    results = utils.nes_data_analysis(forward_bound, reverse_bound, forward_unbound, reverse_unbound)

                    title = self.edgeid_lig_names[edgeid][0] + ' to ' + self.edgeid_lig_names[edgeid][1]

                    # Edge Depiction
                    edge_depiction_string, edge_depiction_image = utils.make_edge_depiction(ligandA, ligandB)

                    # Generate NES floe report
                    report_string = utils.plot_work_pdf(forward_bound,
                                                        reverse_bound,
                                                        forward_unbound,
                                                        reverse_unbound,
                                                        results,
                                                        title,
                                                        edge_depiction_string,
                                                        self.opt['units'])
                    new_record = OERecord()

                    new_record.set_value(Fields.floe_report, report_string)

                    if edge_depiction_string:
                        new_record.set_value(Fields.floe_report_svg_lig_depiction, edge_depiction_string)
                    else:
                        self.opt['Logger'].warn(
                            "It was not possible to generate the edge depiction for the edge id {} "
                            "and it will be skipped".format(edgeid))
                        continue

                    if self.opt['units'] == 'kcal/mol':
                        conv_factor = 4.184
                    else:
                        conv_factor = 1.0

                    label = "BAR score:<br>{:.2f}  &plusmn; {:.2f} {}".format(results['BAR'][0]/conv_factor, results['BAR'][1]/conv_factor, self.opt['units'])
                    new_record.set_value(Fields.floe_report_label, label)
                    new_record.set_value(Fields.floe_report_sort_string, title)
                    new_record.set_value(Fields.FEC.RBFEC.edgeid, edgeid)
                    new_record.set_value(Fields.FEC.RBFEC.edge_name, title)

                    meta_unit = OEFieldMeta().set_option(Meta.Units.Energy.kJ_per_mol)

                    analysis_rec = OERecord()
                    # analysis_rec.set_value(OEField("DDG_BAR", Types.Float, meta=meta_unit), results['BAR'][0])
                    # analysis_rec.set_value(OEField("dDDG_BAR", Types.Float, meta=meta_unit), results['BAR'][1])

                    analysis_rec.set_value(
                        OEField(Fields.FEC.binding_fe.get_name(), Fields.FEC.binding_fe.get_type(), meta=meta_unit),
                        results['BAR'][0])
                    analysis_rec.set_value(
                        OEField(Fields.FEC.binding_fe_err.get_name(), Fields.FEC.binding_fe_err.get_type(),
                                meta=meta_unit), results['BAR'][1])

                    new_record.set_value(Fields.FEC.binding_fe, results['BAR'][0])
                    new_record.set_value(Fields.FEC.binding_fe_err, results['BAR'][1])

                    new_record.set_value(Fields.FEC.RBFEC.NESC.DDG_rec, analysis_rec)

                    work_bound_f = [v for k, v in forward_bound.items()]
                    work_bound_r = [v for k, v in reverse_bound.items()]
                    work_unbound_f = [v for k, v in forward_unbound.items()]
                    work_unbound_r = [v for k, v in reverse_unbound.items()]

                    work_rec = OERecord()
                    work_rec.set_value(OEField("Bound_Forward_Works_OPLMD", Types.FloatVec, meta=meta_unit), work_bound_f)
                    work_rec.set_value(OEField("Bound_Reverse_Works_OPLD", Types.FloatVec,  meta=meta_unit), work_bound_r)
                    work_rec.set_value(OEField("Unbound_Forward_Works_OPLMD", Types.FloatVec, meta=meta_unit), work_unbound_f)
                    work_rec.set_value(OEField("Unbound_Reverse_Works_OPLMD", Types.FloatVec, meta=meta_unit), work_unbound_r)

                    new_record.set_value(Fields.FEC.RBFEC.NESC.work_rec, work_rec)

                    # These Fields are required by the Floe Report and Record Check
                    # to correctly work
                    new_record.set_value(Fields.flask, ligandA)
                    new_record.set_value(Fields.ligand_name, title)
                    new_record.set_value(Fields.title, title + '_' + str(count))
                    new_record.set_value(Fields.ligid, count)
                    new_record.set_value(Fields.confid, 0)

                    # Set the collections
                    new_record.set_value(Fields.collections, self.collections)

                    self.success.emit(new_record)
                    self.opt['Logger'].info("Edge {} Done".format(edgeid), flush=True)
                    count += 1
                except Exception as e:
                    self.opt['Logger'].warn("Error detected skip edge {}".format(edgeid),
                                            str(e), flush=True)
                    continue

        except Exception as e:
            self.opt['Logger'].error("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
        return


class PlotNESResults(RecordPortsMixin, ComputeCube):
    title = "RBFE Plot"
    
    classification = [["RBFE Analysis"]]
    tags = ['Complex', 'Protein', 'Ligand']
    description = """
    This cube generates plots for the Relative binding affinities and the 
    estimated Absolute Binding Affinities versus the provided experimental data. 
    Relevant statistics are calculated and tabled.
    
    Input:
    -------
    Data Record Stream - Streamed-in of records containing predicted relative 
    and absolute binding affinities. 
     
    Output:
    -------
    Floe Report
    """

    uuid = "a62fd733-132b-4619-bb8b-68f373020a79"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    symmetrize = parameters.BooleanParameter(
        'symmetrize',
        default=True,
        help_text="""Select if symmetrize the Relative Binding affinity plot"""
    )

    units = parameters.StringParameter(
        'units',
        choices=['kcal/mol', 'kJ/mol'],
        default='kcal/mol',
        help_text='Units to use to display the plots'
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.edge_pred_dic = dict()
        self.pred_affinity_dic = dict()
        self.exp_affinity_dic = dict()

        if in_orion():
            job_id = environ.get('ORION_JOB_ID')
            self.floe_report = FloeReport.start_report("Affinity Report", job_id=job_id)
        else:
            self.floe_report = LocalFloeReport.start_report("Affinity Report")

    def process(self, record, port):
        try:
            if not record.has_field(Fields.FEC.RBFEC.edge_name):
                raise ValueError("The current record is missing the edge name field")

            edge_name = record.get_value(Fields.FEC.RBFEC.edge_name)
            ligA_name = edge_name.split()[0]
            ligB_name = edge_name.split()[2]

            # Predicted relative binding affinity in kJ/mol
            if not record.has_field(Fields.FEC.RBFEC.NESC.DDG_rec):
                raise ValueError("The current record is missing the Binding Affinity Record")

            DDG_rec = record.get_value(Fields.FEC.RBFEC.NESC.DDG_rec)

            # Free energy values
            if DDG_rec.has_field(Fields.FEC.binding_fe):
                DDG_A_to_B_pred = DDG_rec.get_value(Fields.FEC.binding_fe)
                ddG_A_to_B_pred = DDG_rec.get_value(Fields.FEC.binding_fe_err)
            # NOTE WELL: this elif is only for temporary back compatibility in 1Q2021; remove later
            elif DDG_rec.has_field(OEField('DDG_BAR', Types.Float)):
                DDG_A_to_B_pred = DDG_rec.get_value(OEField('DDG_BAR', Types.Float,
                                                    meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)))
                ddG_A_to_B_pred = DDG_rec.get_value(OEField('dDDG_BAR', Types.Float))
                # End NOTE WELL
            else:
                raise ValueError("The DDG_rec record is missing the Binding Free Energy")

            self.edge_pred_dic[edge_name] = [DDG_A_to_B_pred, ddG_A_to_B_pred]

            if record.has_field(OEField("OPLMD_StateA_ABFE", Types.Record)):

                aff_recA = record.get_value(OEField("OPLMD_StateA_ABFE", Types.Record))

                if ligA_name not in self.pred_affinity_dic:

                    if aff_recA.has_field(Fields.FEC.binding_fe) and aff_recA.has_field(Fields.FEC.binding_fe_err):
                        DG = aff_recA.get_value(Fields.FEC.binding_fe)
                        dG = aff_recA.get_value(Fields.FEC.binding_fe_err)
                        self.pred_affinity_dic[ligA_name] = [DG, dG]

                if ligA_name not in self.exp_affinity_dic:

                    if aff_recA.has_field(Fields.FEC.binding_exptl_fe) and aff_recA.has_field(Fields.FEC.binding_exptl_fe_err):
                        DG_exp = aff_recA.get_value(Fields.FEC.binding_exptl_fe)
                        dG_exp = aff_recA.get_value(Fields.FEC.binding_exptl_fe_err)

                        self.exp_affinity_dic[ligA_name] = [DG_exp, dG_exp]

            if record.has_field(OEField("OPLMD_StateB_ABFE", Types.Record)):

                aff_recB = record.get_value(OEField("OPLMD_StateB_ABFE", Types.Record))

                if ligB_name not in self.pred_affinity_dic:

                    if aff_recB.has_field(Fields.FEC.binding_fe) and aff_recB.has_field(Fields.FEC.binding_fe_err):

                        DG = aff_recB.get_value(Fields.FEC.binding_fe)
                        dG = aff_recB.get_value(Fields.FEC.binding_fe_err)
                        self.pred_affinity_dic[ligB_name] = [DG, dG]

                if ligB_name not in self.exp_affinity_dic:

                    if aff_recB.has_field(Fields.FEC.binding_exptl_fe) and aff_recB.has_field(
                            Fields.FEC.binding_exptl_fe_err):

                        DG_exp = aff_recB.get_value(Fields.FEC.binding_exptl_fe)
                        dG_exp = aff_recB.get_value(Fields.FEC.binding_exptl_fe_err)

                        self.exp_affinity_dic[ligB_name] = [DG_exp, dG_exp]

            # self.success.emit(record)
        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

    def end(self):

        try:
            self.opt['Logger'].info("....Generating Floe Report")

            if self.pred_affinity_dic is None:
                raise ValueError("value_error: no ABFEs were predicted; probably graph weakly connected")

            if not self.exp_affinity_dic:
                self.exp_affinity_dic = None

            # Add expt edge DDG data (where available) to predicted edges
            edge_pred_expt_dic = flrpt.CombinePredExptDDGs(self.exp_affinity_dic, self.edge_pred_dic)

            if self.exp_affinity_dic:

                lig_pred_expt_dic = {name: self.pred_affinity_dic[name] + self.exp_affinity_dic[name] for name in self.pred_affinity_dic if name in self.exp_affinity_dic}

            else:
                lig_pred_expt_dic = {name: self.pred_affinity_dic[name] + [None, None] for name in
                                     self.pred_affinity_dic}

            report_html_str = utils.generate_plots_and_stats(lig_pred_expt_dic,
                                                             edge_pred_expt_dic,
                                                             DDG_symmetrize=self.opt['symmetrize'],
                                                             units=self.opt['units'])

            index = self.floe_report.create_page("index", is_index=True)

            index.set_from_string(report_html_str)

            self.floe_report.finish_report()

        except Exception as e:
            self.opt['Logger'].warn("It was not possible to generate the floe report: {}".format(str(e)))

        return


class PredictDGFromDDG(RecordPortsMixin, ComputeCube):
    title = "PredictDG FromDDG Plot"
    
    classification = [["FEC Analysis"]]
    tags = ['Protein', 'Ligand', 'FEC', 'RBFE', 'NES']
    description = """
    This cube applies the Maximum likelihood Estimator method to predict the absolute
    binding affinities starting from the relative binding affinities (rbfe). If rbfe
    edges are not well graph connected the estimate will not be provided.  
    
    Inputs:
    -------
    Data Record Stream - Streamed-in of records containing predicted relative binding affinities 
    
    Data Record Stream - Streamed-in of records containing ligand information usually provided 
    running the Short Trajectory MD floe
     
    Output:
    -------
    Floe Report
    
    Data record Stream - Streamed-out of records containing per ligand predicted absolute binding 
    affinities. 
    
    Data Record Stream - Streamed-out of records containing  predicted relative 
    and absolute binding affinities used to feed the plot affinity cube
    """

    uuid = "9d8e9390-eb0d-4b97-afec-b0f6482f2843"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    lig_exp_file = FileInputParameter("lig_exp_file", title="OPTIONAL - Ligand Experimental file results",
                                      description="OPTIONAL - Ligand Experimental Results",
                                      required=False,
                                      default=None)

    units = parameters.StringParameter(
        'units',
        choices=['kcal/mol', 'kJ/mol'],
        default='kcal/mol',
        help_text='Units to use for the ABFE output record values'
    )

    graph_port = RecordOutputPort("protein_port")
    bound_port = RecordInputPort("bound_port", initializer=True)

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.exp_dic = dict()
        self.edge_pred_dic = dict()
        self.records = dict()
        self.bound_records = dict()

        for record in self.bound_port:
            if not record.has_field(Fields.ligand_name):
                raise ValueError("The current record is missing the ligand name field")

            lig_name = record.get_value(Fields.ligand_name)
            self.bound_records[lig_name] = record

        file = list(self.args.lig_exp_file)

        if len(file) > 0:
            try:
                for file_obj in file:
                    with TemporaryPath() as path:
                        file_obj.copy_to(path)
                        with open(path, "r") as f:
                            lig_list = f.readlines()
            except Exception as e:
                raise IOError("There was an Error reading the Experimental file: {}".format(str(e)))

            for lig_ln in lig_list:

                if not lig_ln:
                    continue

                if lig_ln.startswith(";") or lig_ln.startswith("#"):
                    continue

                if lig_ln == '\n':
                    continue

                expr = utils.exp_file_grammar(lig_ln)

                lig_name = expr[0]

                if lig_name in self.exp_dic.keys():
                    raise ValueError("Ligand name must be unique: Detected multiple names for: {}".format(lig_name))

                # convert the experimental binding affinity and its error in kJ/mol
                data = []

                if expr[-1] == 'kcal/mol':
                    data.append(expr[1] * 4.184)
                else:
                    data.append(expr[1])

                if len(expr) == 4:
                    if expr[-1] == 'kcal/mol':
                        data.append(expr[2] * 4.184)
                    else:
                        data.append(expr[2])
                else:  # Set the experimental error to zero kJ/mol if not provided
                    data.append(0.0)

                # lig_name_dic['lig_name'] : [DGexp, dGexp]
                self.exp_dic[lig_name] = data

            if self.exp_dic:
                self.opt['Logger'].info("Experimental dictionary built")
            else:
                self.exp_dic = None
                self.opt['Logger'].warn("Experimental dictionary not built. Probably issues with the provided file")

        else:
            self.opt['Logger'].warn("Experimental file not provided")
            self.exp_dic = None

    def process(self, record, port):
        try:
            if not record.has_field(Fields.FEC.RBFEC.edge_name):
                raise ValueError("The current record is missing the edge name field")

            edge_name = record.get_value(Fields.FEC.RBFEC.edge_name)

            # RBFE results subrecord with Predicted relative binding affinity in kJ/mol
            if not record.has_field(Fields.FEC.RBFEC.NESC.DDG_rec):
                raise ValueError("The current record is missing the Binding Affinity Record")

            DDG_rec = record.get_value(Fields.FEC.RBFEC.NESC.DDG_rec)

            # Get calculated relative binding Free energy values
            if DDG_rec.has_field(Fields.FEC.binding_fe):
                DDG_A_to_B_pred = DDG_rec.get_value(Fields.FEC.binding_fe)
                ddG_A_to_B_pred = DDG_rec.get_value(Fields.FEC.binding_fe_err)
            # NOTE WELL: this elif is only for temporary back compatibility in 1Q2021; remove later
            elif DDG_rec.has_field(OEField('DDG_BAR', Types.Float)):
                DDG_A_to_B_pred = DDG_rec.get_value(OEField('DDG_BAR', Types.Float,
                                                            meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)))
                ddG_A_to_B_pred = DDG_rec.get_value(OEField('dDDG_BAR', Types.Float))
                # End NOTE WELL
            else:
                raise ValueError("The DDG_rec record is missing the Binding Free Energy")

            self.edge_pred_dic[edge_name] = [DDG_A_to_B_pred, ddG_A_to_B_pred]
            self.records[edge_name] = record

            # self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

    def end(self):

        # try:
        self.opt['Logger'].info("....Predicting ABFEs from RBFEs")

        # Internally all the RBFE/ABFE results are in kJ/mol
        affinity_dic = utils.predictDGsfromDDGs(self.edge_pred_dic)

        if affinity_dic is None and self.exp_dic is None:
            lig_pred_expt_affinity_dic = {name: [None, None, None, None] for name in self.bound_records}
        else:
            # example lig_pred_expt_affinity_dic[lig_name] = [DPred, dPred, DExp, dExp]
            lig_pred_expt_affinity_dic = CombineAndOffsetPredExptDGs(affinity_dic, self.exp_dic)

        stateA_abfe_field = OEField("OPLMD_StateA_ABFE", Types.Record)
        stateB_abfe_field = OEField("OPLMD_StateB_ABFE", Types.Record)

        meta_unit = OEFieldMeta().set_option(Meta.Units.Energy.kJ_per_mol)

        binding_fe_pred_field = OEField(Fields.FEC.binding_fe.get_name(),
                                        Fields.FEC.binding_fe.get_type(), meta=meta_unit)
        binding_fe_pred_err_field = OEField(Fields.FEC.binding_fe_err.get_name(),
                                            Fields.FEC.binding_fe_err.get_type(), meta=meta_unit)

        binding_fe_exp_field = OEField(Fields.FEC.binding_exptl_fe.get_name(),
                                       Fields.FEC.binding_exptl_fe.get_type(), meta=meta_unit)
        binding_fe_exp_err_field = OEField(Fields.FEC.binding_exptl_fe_err.get_name(),
                                           Fields.FEC.binding_exptl_fe_err.get_type(), meta=meta_unit)

        for lig_name in self.bound_records:

            if lig_name in lig_pred_expt_affinity_dic:

                if self.args.units == 'kcal/mol':
                    meta_unit = OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol)
                    binding_fe_pred_field.set_meta(meta_unit)
                    binding_fe_pred_err_field.set_meta(meta_unit)
                    conv_factor = 4.184
                else:
                    conv_factor = 1.0

                lig_Aff = lig_pred_expt_affinity_dic[lig_name][0]
                lig_Aff_err = lig_pred_expt_affinity_dic[lig_name][1]

                if lig_Aff is not None and lig_Aff_err is not None:
                    bound_rec = self.bound_records[lig_name]
                    bound_rec.set_value(binding_fe_pred_field, lig_pred_expt_affinity_dic[lig_name][0]/conv_factor)
                    bound_rec.set_value(binding_fe_pred_err_field, lig_pred_expt_affinity_dic[lig_name][1]/conv_factor)
                    self.success.emit(bound_rec)

        for edge_name, record in self.records.items():

            # self.opt['Logger'].info("DGG to DG record: {}".format(edge_name))

            ligA_name = edge_name.split()[0]
            ligB_name = edge_name.split()[2]

            if ligA_name in lig_pred_expt_affinity_dic:

                aff_A_rec = OERecord()

                pred_aff_A = lig_pred_expt_affinity_dic[ligA_name][0]
                pred_aff_A_err = lig_pred_expt_affinity_dic[ligA_name][1]

                if pred_aff_A is not None and pred_aff_A_err is not None:

                    aff_A_rec.set_value(binding_fe_pred_field, pred_aff_A)
                    aff_A_rec.set_value(binding_fe_pred_err_field, pred_aff_A_err)

                exp_aff_A = lig_pred_expt_affinity_dic[ligA_name][2]
                exp_aff_A_err = lig_pred_expt_affinity_dic[ligA_name][3]

                if exp_aff_A is not None and exp_aff_A_err is not None:
                    aff_A_rec.set_value(binding_fe_exp_field, exp_aff_A)
                    aff_A_rec.set_value(binding_fe_exp_err_field, exp_aff_A_err)

                if pred_aff_A is not None or exp_aff_A is not None:
                    record.set_value(stateA_abfe_field, aff_A_rec)

            if ligB_name in lig_pred_expt_affinity_dic:

                aff_B_rec = OERecord()

                pred_aff_B = lig_pred_expt_affinity_dic[ligB_name][0]
                pred_aff_B_err = lig_pred_expt_affinity_dic[ligB_name][1]

                if pred_aff_B is not None and pred_aff_B_err is not None:
                    aff_B_rec.set_value(binding_fe_pred_field, pred_aff_B)
                    aff_B_rec.set_value(binding_fe_pred_err_field, pred_aff_B_err)

                exp_aff_B = lig_pred_expt_affinity_dic[ligB_name][2]
                exp_aff_B_err = lig_pred_expt_affinity_dic[ligB_name][3]

                if exp_aff_B is not None and exp_aff_B_err is not None:
                    aff_B_rec.set_value(binding_fe_exp_field, exp_aff_B)
                    aff_B_rec.set_value(binding_fe_exp_err_field, exp_aff_B_err)

                if pred_aff_B is not None or exp_aff_B is not None:
                    record.set_value(stateB_abfe_field, aff_B_rec)

            # else:
            #     self.opt['Logger'].warn("It was not possible to generate the absolute "
            #                             "affinity values for the ligands in the edge: {}".format(edge_name))

            self.graph_port.emit(record)

        # except Exception as e:
        #     self.opt['Logger'].warn("It was not possible to generate the floe report: {}".format(str(e)))

        return


class ParallelNESGMXChimera(ParallelMixin, NESGMXChimera):
    title = "Parallel " + NESGMXChimera.title
    description = "(Parallel) " + NESGMXChimera.description
    uuid = "676baf05-0571-4f14-9f84-d5b5a63729c2"


class ParallelNESGMX(ParallelMixin,  NESGMX):
    title = "Parallel " + NESGMX.title
    description = "(Parallel) " + NESGMX.description
    uuid = "b6640594-8a5a-4e05-89f2-679c5a46691c"
