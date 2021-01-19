from unittest import TestCase
from maze import Zeotype, ImperfectZeotype
from ase import Atoms
from ase.io import read
import os

class TestZeotype(TestCase):

    def test_init_with_no_arguments(self):
        my_zeotype = Zeotype()
        # tests inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # tests empty list attributes
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.silent, False)
        self.assertEqual(my_zeotype.zeolite_type, '')
        self.assertEqual(my_zeotype.site_to_atom_indices, None)
        self.assertEqual(my_zeotype.atom_indices_to_site, None)
        self.assertEqual(my_zeotype.name, 'pristine')
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
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.silent, True)
        self.assertEqual(my_zeotype.zeolite_type, 'friendly')
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
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.silent, False)
        self.assertEqual(my_zeotype.zeolite_type, '')
        self.assertEqual(my_zeotype.site_to_atom_indices, None)
        self.assertEqual(my_zeotype.atom_indices_to_site, None)
        self.assertEqual(my_zeotype.name, 'pristine')
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
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # tests corretly defined parameters
        self.assertEqual(my_zeotype.silent, False)
        self.assertEqual(my_zeotype.zeolite_type, '')
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
        iz = z.get_imperfect_zeolite()

        # tests inheritance
        self.assertIsInstance(iz, ImperfectZeotype)
        self.assertIsInstance(iz, Zeotype)
        self.assertIsInstance(iz, Atoms)
        # tests empty list attributes
        self.assertIsEmptyList(iz.adsorbates)
        self.assertIsEmptyList(iz.sites)
        self.assertIsEmptyList(iz.clusters)
        # tests corretly defined parameters
        self.assertEqual(iz.silent, False)
        self.assertEqual(iz.zeolite_type, '')
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

    def test_count_atomtypes(self):
        self.fail()

    def test_add_cluster(self):
        self.fail()

    def test_integrate_cluster(self):
        self.fail()

    def assertIsEmptyList(self, my_list):
        """
        Args:
            my_list:
        """
        self.assertIsInstance(my_list, list)
        self.assertEqual(len(my_list), 0)
