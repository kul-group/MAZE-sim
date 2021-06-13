from unittest import TestCase
from maze.extra_framework_maker import ExtraFrameworkMaker
from maze.io_zeolite import save_zeolites

# in progress
# class TestOpenDefect(TestCase):
#     def test_OpenDefect(self):
#         input_filepath = 'zeolite_output/open_defect'
#         efm = ExtraFrameworkMaker('BEA')
#         efm.make_extra_frameworks()
#         zeolite_list = []
#         zeolite_list.extend(efm.traj_1Al)
#         zeolite_list.extend(efm.traj_2Al)
#         zeolite_list.append(efm.traj_1Al[0].parent_zeotype)
#         save_zeolites(input_filepath, zeolite_list, zip=False)