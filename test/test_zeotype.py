from unittest import TestCase
from source.zeotype import Zeotype
from ase import Atoms
from ase.io import read

class TestZeotype(TestCase):

    def test_init_with_no_arguments(self):
        my_zeotype = Zeotype()
        # test inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # test empty list attributes
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # test corretly defined parameters
        self.assertEqual(my_zeotype.silent, False)
        self.assertEqual(my_zeotype.zeolite_type, '')
        self.assertEqual(my_zeotype.site_to_atom_indices, None)
        self.assertEqual(my_zeotype.atom_indices_to_site, None)
        self.assertEqual(my_zeotype.name, 'pristine')
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)

    def test_init_with_all_arguments(self):
        my_zeotype = Zeotype(symbols='H', positions=[[0, 0, 10]], numbers=None, tags=[3], momenta=None,
                             masses=None, magmoms=None, charges=None, scaled_positions=None, cell=None,
                             pbc=None, celldisp=None, constraint=None, calculator=None, info=None,
                             velocities=None, silent=True, zeolite_type='friendly',
                             site_to_atom_indices={'T1': 0}, atom_indices_to_site={0: 'T1'}, name='good_friend')
        # TODO: Add in reasonable numbers, tags, momenta, ect.
        # test inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # test empty list attributes
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # test corretly defined parameters
        self.assertEqual(my_zeotype.silent, True)
        self.assertEqual(my_zeotype.zeolite_type, 'friendly')
        self.assertEqual(my_zeotype.site_to_atom_indices, {'T1': 0})
        self.assertEqual(my_zeotype.atom_indices_to_site, {0: 'T1'})
        self.assertEqual(my_zeotype.name, 'good_friend')
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)

    def test_init_with_atoms_obj(self):
        my_atoms = read('GOO.cif')
        my_zeotype = Zeotype(my_atoms)
        # test inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # test empty list attributes
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # test corretly defined parameters
        self.assertEqual(my_zeotype.silent, False)
        self.assertEqual(my_zeotype.zeolite_type, '')
        self.assertEqual(my_zeotype.site_to_atom_indices, None)
        self.assertEqual(my_zeotype.atom_indices_to_site, None)
        self.assertEqual(my_zeotype.name, 'pristine')
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        #test atoms are there and behave the same
        self.assertCountEqual(my_zeotype.get_tags(), my_atoms.get_tags())
        self.assertCountEqual(my_zeotype.get_chemical_symbols(), my_atoms.get_chemical_symbols())

    def test_init_with_zeolite_obj(self):
        z = Zeotype.read_cif_note_sites('G00.cif')
        my_zeotype = Zeotype(z)
        # test inheritance
        self.assertIsInstance(my_zeotype, Zeotype)
        self.assertIsInstance(my_zeotype, Atoms)
        # test empty list attributes
        self.assertIsEmptyList(my_zeotype.adsorbates)
        self.assertIsEmptyList(my_zeotype.sites)
        self.assertIsEmptyList(my_zeotype.clusters)
        # test corretly defined parameters
        self.assertEqual(my_zeotype.silent, False)
        self.assertEqual(my_zeotype.zeolite_type, '')
        self.assertEqual(my_zeotype.site_to_atom_indices, z.site_to_atom_indices)
        self.assertEqual(my_zeotype.atom_indices_to_site, z.atom_indices_to_site)
        self.assertEqual(my_zeotype.name, z.name)
        self.assertEqual(my_zeotype.parent_zeotype, my_zeotype)
        #test atoms are there and behave the same
        self.assertCountEqual(my_zeotype.get_tags(), z.get_tags())
        self.assertCountEqual(my_zeotype.get_chemical_symbols(), z.get_chemical_symbols())

    def test_count_elements(self):
        self.fail()

    def test_count_atomtypes(self):
        self.fail()

    def test_add_cluster(self):
        self.fail()

    def test_integrate_cluster(self):
        self.fail()

    def assertIsEmptyList(self, my_list):
        self.assertIsInstance(my_list, list)
        self.assertEqual(len(my_list), 0)
