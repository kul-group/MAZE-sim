from unittest import TestCase
from maze.zeotypes import ImperfectZeotype, Zeotype
from maze.index_mapper import IndexMapper
from ase import Atoms
import numpy as np
from ase.visualize import view


class TestImperfectZeotype(TestCase):
    def test_build_cap_atoms(self):
        with self.subTest(msg="testing simple case"):
            cap_atoms_dict = {'H': [np.array([0, 0, 0])],
                              'O': [np.array([1, 2, 3])],
                              'N': [np.array([1, 3, 5])]}
            cap_atoms = ImperfectZeotype.build_cap_atoms(cap_atoms_dict)
            for atom in cap_atoms:
                self.assertIn(atom.symbol, cap_atoms_dict)
                self.assertCountEqual(cap_atoms_dict[atom.symbol][0], atom.position)
        with self.subTest(msg="testing multi atoms of same type"):
            cap_atoms_dict = {'H': [np.array([0, 0, 0]), np.array([1, 1, 1])],
                              'O': [np.array([3, 3, 3]), np.array([3, 4, 5])]}
            cap_atoms = ImperfectZeotype.build_cap_atoms(cap_atoms_dict)
            # this complicated for loop checks that the position is in the
            for atom in cap_atoms:
                array_in = False
                for pos in cap_atoms_dict[atom.symbol]:
                    all_equal = True
                    for dict_el, atom_el in zip(pos, atom.position):
                        if dict_el != atom_el:
                            all_equal = False
                            break
                    if all_equal:
                        array_in = True
                        break
                self.assertTrue(array_in)


    def test_cap_atom(self):
        with self.subTest(msg="test hydrogen capping"):
            iz = ImperfectZeotype(Zeotype('O3SiOSi', positions=[[0,0,0], [0, 0, -1], [0, 0, 1]]))
            iz = iz.delete_atoms(2)  # delete Si
            iz = iz.cap_atoms()
            x = 0
            #self.assertEqual(len(iz), 3)
            for atom in iz:
                if atom.symbol == "H":
                    self.assertTrue(np.all(atom.position == np.array([0, 0, 1])))
                if atom.symbol == "O":
                    self.assertTrue(np.all(atom.position == np.array([0, 0, -1])))

    def test_get_cluster(self):
        iz = ImperfectZeotype.make('BEA')
        cluster, od = iz.get_cluster(174)
        with self.subTest(msg='testing length is valid'):
            self.assertEqual(len(cluster) + len(od), len(iz))
        print(cluster.ztype)

    def test_remove_caps(self):
        self.fail()

    def test_integrate_adsorbate(self):
        self.fail()

    def test_remove_adsorbate(self):
        self.fail

    def test_integrate_other_zeotype(self):
        self.fail()

    def test_change_atom_properties(self):
        self.fail()

    def test_build_h_atoms_cap_dict(self):
        self.fail()

    def test_build_o_atoms_cap_dict(self):
        self.fail()

    def test_build_all_atoms_cap_dict(self):
        self.fail()

    def test_get_h_pos(self):
        self.fail()

    def test_find_missing_atom(self):
        self.fail()

    def test_delete_atoms(self):
        self.fail()

    def test_set_attrs_source(self):
        self.fail()

    def test__change_atoms(self):
        self.fail()

    def test_add_atoms(self):
        self.fail()

    def test_remove_addition(self):
        self.fail()

    def test_create_silanol_defect(self):
        self.fail()

    def test_needs_cap(self):
        self.fail()

    def test_get_oxygen_cap_pos(self):
        self.fail()

    def test_get_hydrogen_cap_pos_simple(self):
        self.fail()
