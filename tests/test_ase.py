"""
This file involves using the built in ase tests for atoms to check each class that derives from atoms
"""

from unittest import TestCase
from maze.zeotypes import Zeotype, ImperfectZeotype, Cluster, OpenDefect
import os
from pathlib import Path
from ase import Atoms
import numpy as np

class test_atoms_subclasses(TestCase):
    atoms_subclass_list = [Atoms, Zeotype, ImperfectZeotype, Cluster, OpenDefect]
    class_names = ['Atoms', 'Zeotype', 'ImperfectZeotype', 'Cluster', 'OpenDefect']
    def test_atoms_angles(self):
        for Atoms in self.atoms_subclass_list:
            with self.subTest('testing %s'%str(type(Atoms))):
                atoms = Atoms(['O', 'H', 'H'], positions=[[0., 0., 0.119262],
                                                          [0., 0.763239, -0.477047],
                                                          [0., -0.763239, -0.477047]])

                # Angle no pbc
                self.assertLess(abs(atoms.get_angle(1, 0, 2) - 104), 1e-3)

                atoms.set_cell([2, 2, 2])

                # Across different pbcs
                atoms.set_pbc([True, False, True])
                atoms.wrap()
                self.assertLess(abs(atoms.get_angle(1, 0, 2, mic=True) - 104), 1e-3)

                # Across all True pbc
                atoms.set_pbc(True)
                atoms.wrap()
                self.assertLess(abs(atoms.get_angle(1, 0, 2, mic=True) - 104), 1e-3)

                # Change Angle
                old = atoms.get_angle(1, 0, 2, mic=False)
                atoms.set_angle(1, 0, 2, -10, indices=[2], add=True)
                new = atoms.get_angle(1, 0, 2, mic=False)
                diff = old - new - 10
                self.assertLess(abs(diff), 10e-3)

                # don't actually change angle using indices
                old = atoms.get_angle(1, 0, 2, mic=False)
                atoms.set_angle(1, 0, 2, -10, indices=[2, 1], add=True)
                new = atoms.get_angle(1, 0, 2, mic=False)
                diff = old - new
                self.assertLess(abs(diff), 10e-3)

                # Simple tetrahedron
                tetra_pos = np.array([[0, 0, 0], [1, 0, 0], [.5, np.sqrt(3) * .5, 0],
                                      [.5, np.sqrt(1 / 3.) * .5, np.sqrt(2 / 3.)]])
                atoms = Atoms(['H', 'H', 'H', 'H'],
                              positions=tetra_pos - np.array([.2, 0, 0]))
                angle = 70.5287793655
                self.assertLess(abs(atoms.get_dihedral(0, 1, 2, 3) - angle), 1e-3)

                atoms.set_cell([3, 3, 3])
                atoms.set_pbc(True)
                atoms.wrap()
                self.assertLess(abs(atoms.get_dihedral(0, 1, 2, 3, mic=True) - angle), 1e-3)

    def test_atoms_distance(self):
        for Atoms in self.atoms_subclass_list:
            with self.subTest("testing atoms subclass %s" % str((type(Atoms)))):
                # Setup a chain of H,O,C
                # H-O Dist = 2
                # O-C Dist = 3
                # C-H Dist = 5 with mic=False
                # C-H Dist = 4 with mic=True
                a = Atoms('HOC', positions=[(1, 1, 1), (3, 1, 1), (6, 1, 1)])
                a.set_cell((9, 2, 2))
                a.set_pbc((True, False, False))

                # Calculate indiviually with mic=True
                self.assertEqual(a.get_distance(0, 1, mic=True), 2)
                self.assertEqual(a.get_distance(1, 2, mic=True), 3)
                self.assertEqual(a.get_distance(0, 2, mic=True), 4)

                # Calculate indiviually with mic=False
                self.assertEqual(a.get_distance(0, 1, mic=False), 2)
                self.assertEqual(a.get_distance(1, 2, mic=False), 3)
                self.assertEqual(a.get_distance(0, 2, mic=False), 5)

                # Calculate in groups with mic=True
                self.assertTrue((a.get_distances(0, [1, 2], mic=True) == [2, 4]).all())

                # Calculate in groups with mic=False
                self.assertTrue((a.get_distances(0, [1, 2], mic=False) == [2, 5]).all())

                # Calculate all with mic=True
                self.assertTrue((a.get_all_distances(mic=True) == [[0, 2, 4],
                                                          [2, 0, 3],
                                                          [4, 3, 0]]).all())

                # Calculate all with mic=False
                self.assertTrue((a.get_all_distances(mic=False) == [[0, 2, 5],
                                                           [2, 0, 3],
                                                           [5, 3, 0]]).all())

                # Scale Distance
                old = a.get_distance(0, 1)
                a.set_distance(0, 1, 0.9, add=True, factor=True)
                new = a.get_distance(0, 1)
                diff = new - 0.9 * old
                self.assertLess(abs(diff), 10e-6)

                # Change Distance
                old = a.get_distance(0, 1)
                a.set_distance(0, 1, 0.9, add=True)
                new = a.get_distance(0, 1)
                diff = new - old - 0.9
                self.assertLess(abs(diff), 10e-6)

    def test_atoms_get_duplicates(self):
        from ase.geometry import get_duplicate_atoms
        for Atoms in self.atoms_subclass_list:
            with self.subTest('Atoms in %s'%str(type(Atoms))):

                at = Atoms('H5', positions=[[0., 0., 0.],
                                            [1., 0., 0.],
                                            [1.01, 0, 0],
                                            [3, 2.2, 5.2],
                                            [0.1, -0.01, 0.1]])

                dups = get_duplicate_atoms(at)
                self.assertTrue(all((dups == [[1, 2]]).tolist()))

                dups = get_duplicate_atoms(at, cutoff=0.2)
                self.assertTrue(all((dups == [[0, 4], [1, 2]]).tolist()))

                get_duplicate_atoms(at, delete=True)
                self.assertEqual(len(at), 4)

                at = Atoms('H3', positions=[[0., 0., 0.],
                                            [1., 0., 0.],
                                            [3, 2.2, 5.2]])

                # test if it works if no duplicates are detected.
                get_duplicate_atoms(at, delete=True)
                dups = get_duplicate_atoms(at)

                self.assertEqual(dups.size, 0)

    def test_atoms_get_position(self):
        for Atoms in self.atoms_subclass_list:
            with self.subTest('Atoms in %s'%str(type(Atoms))):
                pbc = [1, 1, 0]
                cell = [[1, 0, 0], [0, 1, 0], [0, 0, 4]]

                positions = [[-0.1, 1.01, -0.5]]
                positions_wrapped = [[0.9, 0.01, -0.5]]

                atoms = Atoms("H", positions=positions, cell=cell, pbc=pbc)

                def test_positions(atoms=atoms):
                    self.assertTrue(np.allclose(positions, atoms.get_positions()))

                def test_positions_wrapped(atoms=atoms):
                    self.assertTrue(np.allclose(positions_wrapped, atoms.get_positions(wrap=True)))

                def test_wrapped_positions(atoms=atoms):
                    atoms.wrap()
                    self.assertTrue(np.allclose(positions_wrapped, atoms.get_positions()))

                test_positions()
                test_positions_wrapped()
                test_wrapped_positions()

    def test_atoms_getitem(self):
        from warnings import warn
        for Atoms in self.atoms_subclass_list:
            with self.subTest('Atoms in %s' % str(type(Atoms))):
                w = Atoms('H2O',
                          positions=[[2.264, 0.639, 0.876],
                                     [0.792, 0.955, 0.608],
                                     [1.347, 0.487, 1.234]],
                          cell=[3, 3, 3],
                          pbc=True)

                try:
                    print(w[True, False])
                except IndexError:
                    pass
                else:
                    # python 3.4 tests skip warnings
                    # other python tests are strict
                    # warnings will be errors
                    warn('')

                self.assertEqual(w[0, 1], w[True, True, False])
                self.assertEqual(w[0, 1], w[0:2])

    def test_atoms_info_copy(self):
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
                at1 = Atoms('H2', positions=[[0., 0., 0.],
                                             [1., 0., 0.]])

                at1.info['str'] = "str"
                at1.info['int'] = 42

                at2 = Atoms(at1)

                self.assertEqual(at2.info, at1.info)

    def test_atoms_instantiation(self):
        from ase import Atom
        for Atoms in self.atoms_subclass_list:
            with self.subTest('Atoms in %s' % str(type(Atoms))):
                """The documentation says:

                    These three are equivalent:

                    >>> d = 1.104  # N2 bondlength
                    >>> a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
                    >>> a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
                    >>> a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d))])

                so let's check"""

                numbers = [7, 7]
                symbols = ["N", "N"]
                dummy_array = 2 * [3 * [0.0]]

                d = 1.104  # N2 bondlength
                a1 = Atoms("N2", [(0, 0, 0), (0, 0, d)])
                a2 = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
                a3 = Atoms([Atom("N", (0, 0, 0)), Atom("N", (0, 0, d))])

                def test_atoms(atoms1=a1, atoms2=a2, atoms3=a3):
                    self.assertEqual(atoms1, atoms2)
                    self.assertEqual(atoms2, atoms3)

                # test redundant keywords
                def test_symbols(numbers=numbers, symbols=symbols):
                    kw = {"numbers": numbers, "symbols": symbols}
                    _test_keywords(**kw)

                def test_momenta(numbers=numbers, momenta=dummy_array):
                    kw = {"momenta": momenta, "velocities": momenta}
                    _test_keywords(numbers=numbers, **kw)

                def test_positions(numbers=numbers, positions=dummy_array):
                    kw = {"positions": positions, "scaled_positions": positions}
                    _test_keywords(numbers=numbers, **kw)

                def _test_keywords(**kw):
                    was_raised = False
                    try:
                        Atoms(**kw)
                    except Exception as inst:
                        assert isinstance(inst, TypeError), inst
                        was_raised = True

                    assert was_raised

                test_atoms()
                test_symbols()
                test_momenta()
                test_positions()

    def test_atoms(self):
        for Atoms in self.atoms_subclass_list:
            with self.subTest('Atoms in %s'%str(type(Atoms))):
                print(Atoms())
                print(Atoms('H2O'))

