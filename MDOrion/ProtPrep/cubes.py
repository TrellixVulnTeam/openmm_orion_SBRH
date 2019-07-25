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

from MDOrion.ProtPrep.utils import ss_bond_fix


class ProteinSetting(RecordPortsMixin, ComputeCube):
    title = "Protein Setting"
    version = "0.1.0"
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

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count=0

    def process(self, record, port):
        try:

            if self.count > 0 and not self.opt['multiple_protein']:
                raise ValueError("Multiple Proteins have been Detected")

            if not record.has_value(Fields.primary_molecule):
                raise ValueError("Missing Primary Molecule field")

            protein = record.get_value(Fields.primary_molecule)

            name = self.opt['protein_title']
            if not name:
                titleFirst12 = protein.GetTitle()[0:12]
                if titleFirst12:
                    name = titleFirst12
                else:
                    name = 'protein'

            protein_ss_fix = ss_bond_fix(protein)

            record.set_value(Fields.title, name)
            record.set_value(Fields.wellid, self.count)
            record.set_value(Fields.primary_molecule, protein_ss_fix)
            record.set_value(Fields.well, protein_ss_fix)

            self.count += 1

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
