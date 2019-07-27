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

import traceback

from MDOrion.Standards import Fields

from openeye import oechem

from orionplatform.mixins import RecordPortsMixin

from floe.api import (ParallelMixin,
                      parameter,
                      ComputeCube)

from oeommtools import packmol


from orionclient.types import ShardCollection

from orionclient.session import (in_orion,
                                 APISession)

from os import environ


class IDSettingCube(RecordPortsMixin, ComputeCube):
    title = "Simulation Well ID Setting"
    version = "0.1.0"
    classification = [["Simulation Well Preparation"]]
    tags = ['Simulation', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube sets the integer ID for each simulation well as well as a descriptive
    title string. If the input molecule 
    on a record has multiple conformers these are split into singles each with 
    its own ID. If a complex will be formed, this cube should be used on ligands
    before forming the complex.
    
    Input:
    -------
    Data record Stream - Streamed input of ligands, one per record

    Output:
    -------
    Data record Stream - Streamed output of records, one per conformer, with title and ID.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.total_count = 0
        self.ligid = -1

    def process(self, record, port):
        try:
            if not record.has_value(Fields.well):
                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("Primary Molecule is missing")
                well = record.get_value(Fields.primary_molecule)

                record.set_value(Fields.well, well)

            well = record.get_value(Fields.well)

            # There should be a ligid; if not, increment the last one
            if not record.has_value(Fields.ligid):
                self.ligid += 1
                record.set_value(Fields.ligid, self.ligid)

            if well.NumConfs() > 1:
                self.opt['Logger'].info("[{}] The well {} has multiple conformers. Each single conformer "
                                        "will be treated as a new molecule".format(self.title,
                                                                                   well.GetTitle()))

            name = well.GetTitle()[0:12]
            if not name:
                name = 'SYS'

            num_conf_counter = 0
            for conf in well.GetConfs():

                conf_mol = oechem.OEMol(conf)

                well_title = name

                if well.GetMaxConfIdx() > 1:
                    well_title += '_c' + str(num_conf_counter)

                conf_mol.SetTitle(well_title)

                record.set_value(Fields.wellid, self.total_count)
                record.set_value(Fields.confid, num_conf_counter)
                record.set_value(Fields.title, well_title)
                record.set_value(Fields.well, conf_mol)

                num_conf_counter += 1

                self.total_count += 1
                self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class CollectionSetting(RecordPortsMixin, ComputeCube):
    title = "Collection Setting"
    version = "0.1.0"
    classification = [["System Preparation"]]
    tags = ['System', 'Complex', 'Protein', 'Ligand']
    description = """
    This cube set a record collection state in open or closed for safety by
    using the cube bool parameter open. A True value will open the record
    collection enabling the shard writing and deleting. If on the record
    the collection field is not present one will be created.

    Input:
    -------
    Data record Stream - Streamed-in of systems such as ligands

    Output:
    -------
    Data Record Stream - Streamed-out of records each one with associated IDs
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    open = parameter.BooleanParameter(
        'open',
        default=True,
        help_text='Open or Close a Collection')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.collection = None

    def process(self, record, port):
        try:

            if in_orion():

                session = APISession

                if record.has_value(Fields.collection):

                    if self.collection is None:

                        collection_id = record.get_value(Fields.collection)

                        collection = session.get_resource(ShardCollection, collection_id)

                        self.collection = collection

                        if self.opt['open']:

                            self.collection.open()

                else:
                    if self.collection is None:

                        job_id = environ.get('ORION_JOB_ID')

                        self.collection = ShardCollection.create(session, job_id)

                        job_id = environ.get('ORION_JOB_ID')

                        if job_id:
                            session.tag_resource(self.collection, "Job {}".format(job_id))

                    record.set_value(Fields.collection, self.collection.id)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return

    def end(self):
        if in_orion():
            if not self.opt['open']:
                if self.collection is not None:
                    self.collection.close()


class SolvationCube(RecordPortsMixin, ComputeCube):
    title = "Solvation Packmol"
    version = "0.1.0"
    classification = [["System Preparation"]]
    tags = ['Complex', 'Protein', 'Ligand', 'Solvation']
    description = """
    The solvation cube solvates a given solute input system in a
    selected mixture of solvents. The solvents can be specified by
    comma separated smiles strings of each solvent component or
    selected keywords like tip3p for tip3p water geometry. For each
    component the user needs to specify its molar fractions as well.
    The solution can be neutralized by adding counter-ions. In addition,
    the ionic solution strength can be set adding salt. The cube
    requires a record as input with a solute molecule to solvate
    and produces an output record with the solvated solute.


     Input:
    -------
    Data record Stream - Streamed-in of system solutes to solvate

    Output:
    -------
    Data Record Stream - Streamed-out of records each with the solvated
    solute
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    density = parameter.DecimalParameter(
        'density',
        default=1.0,
        help_text="Solution density in g/ml")

    padding_distance = parameter.DecimalParameter(
        'padding_distance',
        default=8.0,
        help_text="The padding distance between the solute and the box edge in A")

    distance_between_atoms = parameter.DecimalParameter(
        'distance_between_atoms',
        default=2.0,
        help_text="The minimum distance between atoms in A")

    solvents = parameter.StringParameter(
        'solvents',
        default='tip3p',
        help_text='Select solvents. The solvents are specified as comma separated smiles strings'
                  'e.g. [H]O[H], C(Cl)(Cl)Cl, CS(=O)C or special keywords like tip3p')

    molar_fractions = parameter.StringParameter(
        'molar_fractions',
        default='1.0',
        help_text="Molar fractions of each solvent components. The molar fractions are specified"
                  "as comma separated molar fractions strings e.g. 0.5,0.2,0.3")

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text='Output Packmol log')

    geometry = parameter.StringParameter(
        'geometry',
        default='box',
        choices=['box', 'sphere'],
        help_text="Geometry selection: box or sphere. Sphere cannot be used as periodic system "
                  "along with MD simulation")

    close_solvent = parameter.BooleanParameter(
        'close_solvent',
        default=False,
        help_text="If Checked/True solvent molecules will be placed very close to the solute")

    salt = parameter.StringParameter(
        'salt',
        default='[Na+], [Cl-]',
        help_text='Salt type. The salt is specified as list of smiles strings. '
                  'Each smiles string is the salt component dissociated in the '
                  'solution e.g. Na+, Cl-')

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=0.0,
        help_text="Salt concentration in millimolar")

    neutralize_solute = parameter.BooleanParameter(
        'neutralize_solute',
        default=True,
        help_text='Neutralize the solute by adding Na+ and Cl- counter-ions based on'
                  'the solute formal charge')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = self.opt

            if not record.has_value(Fields.well):
                raise ValueError("Missing the Well Molecule Field")

            solute = record.get_value(Fields.well)

            if not record.has_value(Fields.title):
                self.log.warn("Missing Title field")
                solute_title = solute.GetTitle()[0:12]
            else:
                solute_title = record.get_value(Fields.title)

            self.log.info("[{}] solvating well {}".format(self.title, solute_title))

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solute) if dp.GetTag() in
                            ["solvents", "molar_fractions", "density"]}
            if new_args:
                for k in new_args:
                    if k == 'molar_fractions':
                        continue
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solute_title, new_args))
                opt.update(new_args)

            # Solvate the system
            sol_system = packmol.oesolvate(solute, **opt)
            self.log.info("[{}] Solvated simulation well {} yielding {} atoms overall".format(self.title,
                                                                    solute_title, sol_system.NumAtoms()))
            sol_system.SetTitle(solute.GetTitle())

            record.set_value(Fields.well, sol_system)
            record.set_value(Fields.title, solute_title)

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return


class ParallelSolvationCube(ParallelMixin, SolvationCube):
    title = "Parallel " + SolvationCube.title
    description = "(Parallel) " + SolvationCube.description
