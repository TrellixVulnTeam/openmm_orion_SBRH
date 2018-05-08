import traceback

from cuberecord import OERecordComputeCube

from Standards import Fields


class ProteinSetting(OERecordComputeCube):
    title = "Protein Setting"
    version = "0.0.0"
    classification = [["Protein Preparation", "OEChem"]]
    tags = ['OEChem']
    description = """
    This cube is setting the protein to perform MD. A title is generated by
    using the protein name 

    Input:
    -------
    Data record Stream - Streamed-in of just one protein

    Output:
    -------
    Data Record Stream - Emits the MD set protein
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.count=0

    def process(self, record, port):
        try:

            if self.count > 0:
                raise ValueError("Multiple Proteins have been Detected "
                                 "Currently just one Protein is supported as input")

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing '{}' field".format(Fields.primary_molecule.get_name()))
                raise ValueError("Missing Primary Molecule")

            protein = record.get_value(Fields.primary_molecule)

            name = protein.GetTitle()[0:12]
            if not name:
                name = 'PRT'

            record.set_value(Fields.title, name)
            self.count += 1

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)