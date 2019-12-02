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

from floe.api import (parameter,
                      ComputeCube)

from orionplatform.mixins import RecordPortsMixin

# from MDOrion.ProtPrep.utils import ss_bond_fix

from openeye import oechem

from oeommtools import utils

from MDOrion.ForceField import utils as ffutils


class ProteinSetting(RecordPortsMixin, ComputeCube):
    title = "Protein Setting"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['Protein']
    description = """
    This cube is currently used just to check that one
    protein as provided as system input to perform MD. A
    cube parameter can be set to change this behaviour but currently
    multiple protein are not supported by other MD cubes and the user 
    currently have not to set this parameter. A field record 
    title is also generated by using the protein name 

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of a single protein system

    Output:
    -------
    oechem.OEDataRecord - Streamed-out of a single protein
    """

    uuid = "66d29c80-493a-4057-88cc-2b93b1018011"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    multiple_protein = parameter.BooleanParameter(
        'multiple_protein',
        default=False,
        help_text="If Checked/True multiple protein will be allowed")

    protein_title = parameter.StringParameter(
        'protein_title',
        default='',
        help_text='Optional replacement for the protein title'
    )

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default=sorted(ffutils.proteinff)[0],
        choices=sorted(ffutils.proteinff),
        help_text='Force field parameters to be applied to the protein')

    solvent_forcefield = parameter.StringParameter(
        'solvent_forcefield',
        default=sorted(ffutils.solventff)[0],
        help_text='Force field parameters to be applied to the water')

    other_forcefield = parameter.StringParameter(
        'other_forcefield',
        default=sorted(ffutils.otherff)[0],
        choices=sorted(ffutils.otherff),
        help_text='Force field used to parametrize other molecules not recognized by the '
                  'protein force field like excipients')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count=0
        self.opt['CubeTitle'] = self.title

    def process(self, record, port):
        try:

            if self.count > 0 and not self.opt['multiple_protein']:
                raise ValueError("Multiple Proteins have been Detected")

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing Primary Molecule field")

            protein = record.get_value(Fields.primary_molecule)

            # Removing Interaction Hint Container, Style and PDB Data
            oechem.OEDeleteInteractionsHintSerializationData(protein)
            oechem.OEDeleteInteractionsHintSerializationIds(protein)
            oechem.OEClearStyle(protein)
            oechem.OEClearPDBData(protein)

            # Check Protein parametrization
            protein_sp, ligand, water, others = utils.split(protein)

            # Parametrization Checking
            if protein_sp.NumAtoms() > 0:
                ffutils.applyffProtein(protein_sp, self.opt)
            else:
                raise ValueError("The protein does not have any atoms")
            if water.NumAtoms() > 0:
                ffutils.applyffWater(water, self.opt)
            if others.NumAtoms() > 0:
                # Unique prefix name used to output parametrization files
                self.opt['prefix_name'] = 'protein' + '_' + str(0)

                ffutils.applyffExcipients(others, self.opt)

            name = self.opt['protein_title']

            if not name:
                title_first12 = protein.GetTitle()[0:12]
                if title_first12:
                    name = title_first12
                else:
                    name = 'protein'

            # protein_ss_fix = ss_bond_fix(protein)

            record.set_value(Fields.title, name)
            record.set_value(Fields.flaskid, self.count)
            # record.set_value(Fields.primary_molecule, protein_ss_fix)
            # record.set_value(Fields.flask, protein_ss_fix)

            record.set_value(Fields.primary_molecule, protein)
            record.set_value(Fields.flask, protein)

            self.count += 1

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
