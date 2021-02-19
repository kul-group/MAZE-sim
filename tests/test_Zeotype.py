from unittest import TestCase
from maze.zeotypes import Zeotype, ImperfectZeotype
from ase import Atoms
from ase.io import read
import os

class TestZeotype(TestCase):

    # testing init methods
    def test_init_with_no_arguments(self):
        my_zeotype = Zeotype()
        # tests inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.site_to_atom_indices, None)
        self.assertEqual(my_zeotype.atom_indices_to_site, None)
        self.assertEqual('parent', my_zeotype.name)
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        self.assertNotEqual(my_zeotype.index_mapper, None)

    def test_init_with_all_arguments(self):
        my_zeotype = Zeotype(symbols='H', positions=[[0, 0, 10]], numbers=None, tags=[3], momenta=None,
                             masses=None, magmoms=None, charges=None, scaled_positions=None, cell=None,
                             pbc=None, celldisp=None, constraint=None, calculator=None, info=None,
                             velocities=None, silent=True, zeolite_type='friendly',
                             site_to_atom_indices={'T1': 0}, atom_indices_to_site={0: 'T1'}, name='good_friend')
        # TODO: Add in reasonable numbers, tags, momenta, ect.
        # tests inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.site_to_atom_indices, {'T1': 0})
        self.assertEqual(my_zeotype.atom_indices_to_site, {0: 'T1'})
        self.assertEqual(my_zeotype.name, 'good_friend')
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        self.assertNotEqual(my_zeotype.index_mapper, None)

    def test_init_with_atoms_obj(self):
        my_atoms = read('GOO.cif')
        my_zeotype = Zeotype(my_atoms)
        # tests inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.site_to_atom_indices, None)
        self.assertEqual(my_zeotype.atom_indices_to_site, None)
        self.assertEqual('parent', my_zeotype.name)
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        #tests atoms are there and behave the same
        self.assertCountEqual(my_zeotype.get_tags(), my_atoms.get_tags())
        self.assertCountEqual(my_zeotype.get_chemical_symbols(), my_atoms.get_chemical_symbols())

    def test_init_with_zeolite_obj(self):
        z = Zeotype.build_from_cif_with_labels('GOO.cif')
        my_zeotype = Zeotype(z)
        print(os.getcwd())
        # tests inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.site_to_atom_indices, z.site_to_atom_indices)
        self.assertEqual(my_zeotype.atom_indices_to_site, z.atom_indices_to_site)
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
        z = Zeotype.build_from_cif_with_labels('GOO.cif')
        iz = z.get_imperfect_zeotype()

        # tests inheritance
        self.assertIsInstance(iz, ImperfectZeotype)
        self.assertIsInstance(iz, Zeotype)
        self.assertIsInstance(iz, Atoms)
        # tests empty list attributes
        # tests corretly defined parameters
        self.assertEqual(iz.site_to_atom_indices, None)
        self.assertEqual(iz.atom_indices_to_site, None)
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


    def test_build_from_cif_with_labels(self):  # TODO: Write test case for building from CIF file (use rel path)
        abs_path = '/Users/dda/Code/MAZE-sim/data/BEA.cif'
        Zeotype.build_from_cif_with_labels(abs_path)
        self.assertEqual(True, True)

    # other methods

    def test_get_available_symbols(self):  # TODO: Write this test case
        self.assertEqual(True, True)

    def test__read_cif_note_siteJan2021Update(self):
        pass

    def test_read_cif_note_sites(self):
        pass

    def test_get_imperfect_zeolite(self):
        pass

    def test_update_nl(self):
        pass

    def test_get_hetero_atoms(self):
        pass

    def test_get_atom_types(self):
        pass

    def test_count_elements(self):
        pass

    def test_get_cluster(self):
        pass

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




