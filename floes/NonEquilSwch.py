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

from os import path

from floe.api import WorkFloe

from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube

from MDOrion.SubFloes.SubfloeFunctions import nes_gmx_subfloe

from MDOrion.Flask.cubes import CollectionSetting

from MDOrion.Flask.cubes import ParallelRecordSizeCheck


floe_title = 'Non-Equilibrium Switching'
tags_for_floe = ['MD', 'FECalc']
#
tag_str = ''.join(' [{}]'.format(tag) for tag in tags_for_floe)
job = WorkFloe(floe_title, title=floe_title+tag_str)
job.classification = [tags_for_floe]
job.tags = tags_for_floe

job.description = open(path.join(path.dirname(__file__), 'NonEquilSwch_desc.rst'), 'r').read()

job.uuid = "74cd690f-f98a-47e0-bfa4-1858e4080dc3"


# Unbound Reader
ibn = DatasetReaderCube("BoundReader", title="Bound Reader")
ibn.promote_parameter("data_in", promoted_name="bound", title="Bound Input Dataset", description="Bound Input Dataset", order=0)

iun = DatasetReaderCube("UnboundReader", title="Unbound Reader")
iun.promote_parameter("data_in", promoted_name="unbound", title="Unbound Input Dataset", description="Unbound Input Dataset", order=1)

# This cube is necessary for the correct work of collection and shard
coll_open = CollectionSetting("OpenNESCollection", title="OpenNESCollection")
coll_open.set_parameters(open=True)
coll_open.set_parameters(write_new_collection='NES_OPLMD')

# This cube is necessary for the correct working of collections and shards
coll_close = CollectionSetting("CloseCollection", title="Close Collection")
coll_close.set_parameters(open=False)

rec_check = ParallelRecordSizeCheck("Record Check Success", title="Record Size Checking")
rec_check_abfe = ParallelRecordSizeCheck("Record Check Success ABFE", title="Affinity Record Size Checking")

ofs_nes = DatasetWriterCube('ofs', title='NES Out')
ofs_nes.promote_parameter("data_out", promoted_name="out",
                          title="NES Dataset Out",
                          description="NES Dataset Out", order=4)

ofs_abfe = DatasetWriterCube('ofs_abfe', title='Affinity Out')
ofs_abfe.promote_parameter("data_out", promoted_name="abfe",
                           title="Affinity Out",
                           description="Affinity Out", order=5)

fail = DatasetWriterCube('fail', title='NES Failures')
fail.promote_parameter("data_out", promoted_name="fail", title="NES Failures",
                       description="NES Dataset Failures out", order=6)

job.add_cubes(iun, ibn, coll_open, coll_close, rec_check, rec_check_abfe, ofs_nes, ofs_abfe, fail)

nes_subfloe_options = dict()
nes_subfloe_options['edge_map_file'] = 'map'
nes_subfloe_options['n_traj_frames'] = 80
nes_subfloe_options['nes_switch_time_in_ns'] = 0.05

input_port_dic = {'input_open_collection_port': coll_open.success,
                  'input_bound_port': ibn.success}
output_port_dic = {'output_nes_port': coll_close.intake,
                   'output_abfe_port': rec_check_abfe.intake,
                   'output_fail_port': rec_check.fail_in}

nes_gmx_subfloe(job, input_port_dic, output_port_dic, nes_subfloe_options)

iun.success.connect(coll_open.intake)
ibn.success.connect(coll_open.intake)
coll_close.success.connect(rec_check.intake)
rec_check.success.connect(ofs_nes.intake)
rec_check_abfe.success.connect(ofs_abfe.intake)

coll_open.failure.connect(rec_check.fail_in)
coll_close.failure.connect(rec_check.fail_in)
rec_check.failure.connect(fail.intake)
rec_check_abfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
