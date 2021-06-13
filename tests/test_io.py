from unittest import TestCase
from maze.zeolite import Zeolite
from maze.io_zeolite import save_zeolites, read_zeolites
import glob
import os
from ase import Atoms
import json
import numpy as np
from ase.visualize import view


class IOZeolites(TestCase):
    def test_read_zeotypes(self):
        output_filepath = 'zeolite_output/test_zeo'
        cha = Zeolite.make('CHA')
        cha2 = cha.delete_atoms(cha.site_to_atom_indices['T1'])
        water = Atoms("HHeO", positions=[[1,1,1], [0,0,0], [-1, -1, -1]])
        #water = Atoms("HO", positions=[[0,0,0], [-1, -1, -1]])

        cha3 = cha2.add_atoms(water, 'water')
        zeolite_list = [cha, cha2, cha3, cha3.parent_zeotype]
        save_zeolites(output_filepath, zeolite_list, zip=False)
        folder_file_list = glob.glob(output_filepath + '/**')

        with self.subTest('test that folders have been made'):
            for z in zeolite_list:
                self.assertIn(os.path.join(output_filepath, z.name), folder_file_list)

        # import pandas as pd
        # print(pd.DataFrame(cha.index_mapper.main_index).T.to_string())

    # def test_read_zeolites(self):
    #     input_filepath = 'zeolite_output/test_zeo'
    #     zeotype_dict = read_zeolites(input_filepath, str_ext='.traj', zipped=False)
    #     with self.subTest('test same index mapper and same parent'):
    #         zeotypes = list(zeotype_dict.values())
    #         for i in range(1, len(zeotypes)):
    #             self.assertEqual(zeotypes[0].index_mapper, zeotypes[i].index_mapper)
    #             self.assertEqual(zeotypes[0].parent_zeotype, zeotypes[i].parent_zeotype)
    #
    #     with open(os.path.join(input_filepath, 'index_mapper.json')) as f:
    #         original_index_mapper = json.load(f)['main_index']
    #     with self.subTest('test index mapper recreated'):
    #         for row_index, row in original_index_mapper.items():
    #             for name, item_index in row.items():
    #                 if item_index is not None:
    #                     item_index = int(item_index)
    #                 self.assertEqual(item_index, zeotypes[0].index_mapper.main_index[int(row_index)][name],
    #                                  msg=f'testing main index {row_index}, name {name}, item_index {item_index}')
    #
    #     # TODO: Add in glob cmd
    #     #import pandas as pd
    #     #print(pd.DataFrame(zeotypes[0].index_mapper.main_index).T.to_string())
    #
    # def test_tag_zeolite(self):
    #     with self.subTest('test overlapping atoms'):
    #         cha1 = Zeolite.make('CHA')
    #         cha2 = Zeolite.make('CHA')
    #         for atom in cha2:
    #             cha2[atom.index].tag = atom.index + 5
    #         self.assertCountEqual(cha1.get_tags(), cha2.get_tags())
    #         for a1, a2 in zip(cha1, cha2):
    #             self.assertEqual(a1.tag, a2.tag)
    #
    #     with self.subTest('test offset atoms'):
    #         cha1 = Zeolite.make('CHA')
    #         cha2 = Zeolite.make('CHA')
    #         for atom in cha2:
    #             cha2[atom.index].position = cha2[atom.index].position + np.random.random(3)
    #             cha2[atom.index].tag = atom.index + 5
    #         self.assertCountEqual(cha1.get_tags(), cha2.get_tags())
    #         for a1, a2 in zip(cha1, cha2):
    #             self.assertEqual(a1.tag, a2.tag)
