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

from floe.api import (parameters,
                      ComputeCube)

from orionplatform.mixins import RecordPortsMixin

from openeye import oechem

from MDOrion.ForceField.utils import MDComponents


class MDSetting(RecordPortsMixin, ComputeCube):
    title = "MD Setting"
    # version = "0.1.4"
    classification = [["System Preparation"]]
    tags = ['Protein']
    description = """
    This cube is used to componentize the starting extended protein.
    The cube detects if a DU is present on the record and will try 
    to extract the components saving them in ad-hoc container. If the 
    DU is not found, the cube will try to create a DU and if it fails 
    the primary molecule present on the record will be componentize. 
    """

    uuid = "b85d652f-188a-4cc0-aefd-35c98e737f8d"

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 14000},
        "spot_policy": {"default": "Prohibited"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    flask_title = parameters.StringParameter(
        'flask_title',
        default='',
        help_text='Flask Title')

    multiple_flasks = parameters.BooleanParameter(
        'multiple_flasks',
        default=False,
        help_text="If Checked/True multiple flasks will be allowed")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count=0
        self.opt['CubeTitle'] = self.title

    def process(self, record, port):
        try:

            if self.count > 0 and not self.opt['multiple_flasks']:
                raise ValueError("Multiple Flasks have been Detected")

            name = self.opt['flask_title']

            if not name:
                name = 'protein'

            if record.has_value(Fields.design_unit_from_spruce):

                du = oechem.OEDesignUnit()

                if not oechem.OEReadDesignUnitFromBytes(du, record.get_value(Fields.design_unit_from_spruce)):
                    raise ValueError("It was not possible Reading the Design Unit from the record")

                self.opt['Logger'].info("[{}] Design Unit Detected".format(self.title))

                md_components = MDComponents(du, components_title=name)

            else:  # The extended protein is already prepared to MD standard

                if not record.has_value(Fields.primary_molecule):
                    raise ValueError("Missing Primary Molecule field")

                molecules = record.get_value(Fields.primary_molecule)

                md_components = MDComponents(molecules, components_title=name)

            self.opt['Logger'].info(md_components.get_info)

            record.set_value(Fields.md_components, md_components)
            record.set_value(Fields.title, name)
            record.set_value(Fields.flaskid, self.count)

            self.count += 1

            self.success.emit(record)

        except Exception as e:

            print("Failed to complete", str(e), flush=True)
            self.opt['Logger'].info('Exception {} {}'.format(str(e), self.title))
            self.log.error(traceback.format_exc())
            self.failure.emit(record)

        return
