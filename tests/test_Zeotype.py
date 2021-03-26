from unittest import TestCase
from maze.zeolite import PerfectZeolite, Zeolite, Cluster, OpenDefect
import ase
import ase.data
from ase import Atoms
from ase.io import read
import os

class TestZeotype(TestCase):

    # testing init methods
    def test_init_with_no_arguments(self):
        my_zeotype = PerfectZeolite()
        # tests inheritance
        self.assertIsInstance(my_zeotype, PerfectZeolite)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype._site_to_atom_indices, None)
        self.assertEqual(my_zeotype._atom_indices_to_site, None)
        self.assertEqual('parent', my_zeotype.name)
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        self.assertNotEqual(my_zeotype.index_mapper, None)

    def test_init_with_all_arguments(self):
        my_zeotype = PerfectZeolite(symbols='H', positions=[[0, 0, 10]], numbers=None, tags=[3], momenta=None,
                                    masses=None, magmoms=None, charges=None, scaled_positions=None, cell=None,
                                    pbc=None, celldisp=None, constraint=None, calculator=None, info=None,
                                    velocities=None, silent=True, zeolite_type='friendly',
                                    site_to_atom_indices={'T1': 0}, atom_indices_to_site={0: 'T1'}, ztype='good_friend')
        # TODO: Add in reasonable numbers, tags, momenta, ect.
        # tests inheritance
        self.assertIsInstance(my_zeotype, PerfectZeolite)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype._site_to_atom_indices, {'T1': 0})
        self.assertEqual(my_zeotype._atom_indices_to_site, {0: 'T1'})
        self.assertEqual(my_zeotype.name, 'good_friend')
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        self.assertNotEqual(my_zeotype.index_mapper, None)

    def test_init_with_atoms_obj(self):
        my_atoms = read('GOO.cif')
        my_zeotype = PerfectZeolite(my_atoms)
        # tests inheritance
        self.assertIsInstance(my_zeotype, PerfectZeolite)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype._site_to_atom_indices, None)
        self.assertEqual(my_zeotype._atom_indices_to_site, None)
        self.assertEqual('parent', my_zeotype.name)
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        #tests atoms are there and behave the same
        self.assertCountEqual(my_zeotype.get_tags(), my_atoms.get_tags())
        self.assertCountEqual(my_zeotype.get_chemical_symbols(), my_atoms.get_chemical_symbols())

    def test_init_with_zeolite_obj(self):
        z = PerfectZeolite.build_from_cif_with_labels('GOO.cif')
        my_zeotype = PerfectZeolite(z)
        print(os.getcwd())
        # tests inheritance
        self.assertIsInstance(my_zeotype, PerfectZeolite)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype._site_to_atom_indices, z._site_to_atom_indices)
        self.assertEqual(my_zeotype._atom_indices_to_site, z._atom_indices_to_site)
        self.assertEqual(my_zeotype.name, z.name)
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        #tests atoms are there and behave the same
        self.assertCountEqual(my_zeotype.get_tags(), z.get_tags())
        self.assertCountEqual(my_zeotype.get_chemical_symbols(), z.get_chemical_symbols())
        self.assertNotEqual(id(my_zeotype.parent_zeotype), id(z.parent_zeotype))
        self.assertNotEqual(id(my_zeotype.index_mapper), id(z.index_mapper))
        self.assertNotEqual(id(my_zeotype.index_mapper), id(z.parent_zeotype))
        self.assertNotEqual(my_zeotype.index_mapper, None)

    def test_build_iz_from_z(self):
        z = PerfectZeolite.build_from_cif_with_labels('GOO.cif')
        iz = Zeolite(z)
        # tests inheritance
        self.assertIsInstance(iz, Zeolite)
        self.assertIsInstance(iz, PerfectZeolite)
        self.assertIsInstance(iz, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(iz._site_to_atom_indices, None)
        self.assertEqual(iz._atom_indices_to_site, None)
        #self.assertEqual(iz.name)
        self.assertEqual(iz.parent_zeotype, z)
        #tests atoms are there and behave the same
        self.assertCountEqual(iz.get_tags(), z.get_tags())
        self.assertCountEqual(iz.get_chemical_symbols(), z.get_chemical_symbols())
        self.assertEqual(id(iz.parent_zeotype), id(z.parent_zeotype))
        self.assertEqual(id(iz.index_mapper), id(z.index_mapper))
        self.assertEqual(id(iz.parent_zeotype), id(z.parent_zeotype))
        self.assertNotEqual(iz.index_mapper, None)

        self.assertEqual(iz.parent_zeotype, z)
        self.assertEqual(iz.index_mapper, z.index_mapper)
        self.assertEqual(iz.index_mapper, z.index_mapper)


    def test_make(self):
        z = PerfectZeolite.make('CHA')

        cha_site_dict = {'O1': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                         'O2': [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
                         'O3': [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53],
                         'O4': [54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71],
                         'T1': [72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                                91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107]}

        self.assertIsInstance(z, PerfectZeolite)
        self.assertIs(z, z.parent_zeotype)
        self.assertDictEqual(cha_site_dict, z._site_to_atom_indices)

    def test_build_from_cif_with_labels(self):  # TODO: Write test case for building from CIF file (use rel path)
        abs_path = '/Users/dda/Code/MAZE-sim/tests/data/CHA.cif'
        cha_site_dict = {'O1': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                         'O2': [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
                         'O3': [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53],
                         'O4': [54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71],
                         'T1': [72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                                91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107]}
        z = PerfectZeolite.build_from_cif_with_labels(abs_path)
        self.assertIsInstance(z, PerfectZeolite)
        self.assertIs(z, z.parent_zeotype)
        self.assertDictEqual(cha_site_dict, z._site_to_atom_indices)

    # other methods

    def test_get_available_symbols(self):  # TODO: Write this test case
        all_symbols = ase.data.chemical_symbols
        atoms_to_exclude = ['H', 'O', 'Ca', 'He']
        with self.subTest(msg='All but these four'):
            available_syms = PerfectZeolite._get_available_symbols(atoms_to_exclude)
            self.assertEqual(len(available_syms), len(all_symbols) - len(atoms_to_exclude))
            for sym in atoms_to_exclude:
                self.assertNotIn(sym, available_syms)
                self.assertIn(sym, all_symbols)


    def test__read_cif_note_siteJan2021Update(self):
        pass

    def test_read_cif_note_sites(self):
        pass

    def test_get_imperfect_zeolite(self):
        z = PerfectZeolite.make('CHA')
        iz = z.get_imperfect_zeotype()
        self.assertIsInstance(iz, Zeolite)
        self.assertIs(iz.parent_zeotype, z)
        self.assertIs(iz.index_mapper, z.index_mapper)

    def test_update_nl(self):
        pass

    def test_get_hetero_atoms(self):
        z = PerfectZeolite.make('CHA')
        with self.subTest(msg='testing default arguments'):
            z[3].symbol = 'Hf'
            z[5].symbol = 'Sn'
            self.assertCountEqual([3,5],z.get_hetero_atoms())
        with self.subTest(msg='Custom hetero atoms list'):
            # first get a list of O atoms
            O_list = []
            for atom in z:
                if atom.symbol == 'O':
                    O_list.append(atom.index)
            self.assertCountEqual(O_list, z.get_hetero_atoms(['O']))

    def test_get_atom_types(self):
        z = PerfectZeolite.make('CHA')
        get_atom_types_dict = {'framework-O': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                               20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
                                               38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                                               56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71],
                               'framework-Si': [72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87,
                                                88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102,
                                                103, 104, 105, 106, 107]}

        with self.subTest(msg='test get atom types'):
            self.assertDictEqual(z.get_atom_types(), get_atom_types_dict)


    def test_count_elements(self):
        z = PerfectZeolite('H3O5F15')
        indices, count = z.count_elements()
        indices_dict = {'H': [i for i in range(0, 3)],
                        'O': [i for i in range(3, 3+5)],
                        'F': [i for i in range(8, 8+15)]}

        count_dict = {'H': 3, 'O': 5, 'F': 15}
        self.assertCountEqual(indices_dict.keys(), count_dict.keys())
        for key in indices.keys():
            self.assertCountEqual(indices[key], indices_dict[key])
        self.assertDictEqual(count_dict, count)


    def test_get_cluster(self):
        z = PerfectZeolite.make('CHA')
        cluster, opendefect = z.get_cluster(0, 100, 3)
        self.assertIsInstance(cluster, Cluster)
        self.assertIsInstance(opendefect, OpenDefect)


    def test_find_silanol_groups(self):
        pass

    def test_find_silanol_nest_T_sites(self):
        pass

    def test_get_indices_compliment(self):
        pass

    def test_get_indices(self):
        pass

    def test_get_site_type(self):
        pass

    def test__get_old_to_new_map(self):
        pass

    def test_copy(self):
        pass

    def test__del__(self):
        pass




