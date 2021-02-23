"""
This file involves using the built in ase tests for atoms to check each class that derives from atoms
"""

from unittest import TestCase
from maze.zeotypes import Zeotype, ImperfectZeotype, Cluster, OpenDefect
from maze.adsorbate import Adsorbate
import os
from pathlib import Path
from ase import Atoms
import numpy as np

class test_atoms_subclasses(TestCase):
    atoms_subclass_list = [Atoms, Zeotype, ImperfectZeotype, Cluster, OpenDefect, Adsorbate]
    class_names = ['Atoms', 'Zeotype', 'ImperfectZeotype', 'Cluster', 'OpenDefect', 'Adsorbate']
    def test_atoms_angles(self):
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
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
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
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
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):

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
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
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
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
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
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
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
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
                print(Atoms())
                print(Atoms('H2O'))

    def test_diffusion_coefficient(self):
        from ase.md.analysis import DiffusionCoefficient
        from ase.units import fs as fs_conversion
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):

                eps = 1e-10
                # Creating simple trajectories
                # Textbook case. The displacement coefficient should be 0.5 A^2 / fs except for the final molecule

                ###### He atom

                he = Atoms('He', positions=[(0, 0, 0)])
                traj_he = [he.copy() for i in range(2)]
                traj_he[1].set_positions([(1, 1, 1)])

                timestep = 1 * fs_conversion  # fs

                dc_he = DiffusionCoefficient(traj_he, timestep)
                dc_he.calculate(ignore_n_images=0, number_of_segments=1)
                ans = dc_he.get_diffusion_coefficients()[0][0]
                # Answer in \AA^2/<ASE time unit>
                ans_orig = 5.0e-01 / fs_conversion
                # dc_he.print_data()

                self.assertLess(abs(ans - ans_orig), eps)

                ###### CO molecule

                co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1)])
                traj_co = [co.copy() for i in range(2)]
                traj_co[1].set_positions([(-1, -1, -1), (-1, -1, 0)])

                dc_co = DiffusionCoefficient(traj_co, timestep, molecule=False)
                dc_co.calculate(ignore_n_images=0, number_of_segments=1)
                ans = dc_co.get_diffusion_coefficients()[0][0]
                # dc_co.print_data()

                self.assertLess(abs(ans - ans_orig),  eps)

                for index in range(2):
                    dc_co = DiffusionCoefficient(traj_co, timestep, atom_indices=[index], molecule=False)
                    dc_co.calculate()
                    ans = dc_co.get_diffusion_coefficients()[0][0]
                    assert (abs(ans - ans_orig) < eps)

                dc_co = DiffusionCoefficient(traj_co, timestep, molecule=True)
                dc_co.calculate(ignore_n_images=0, number_of_segments=1)
                ans = dc_co.get_diffusion_coefficients()[0][0]
                # dc_co.print_data()

                self.assertLess(abs(ans - ans_orig), eps)

    def test_fix_bond_length_mic(self):
        from ase.calculators.lj import LennardJones
        from ase.constraints import FixBondLength
        from ase.optimize import FIRE
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):

                for wrap in [False, True]:
                    a = Atoms('CCC',
                                  positions=[[1, 0, 5],
                                             [0, 1, 5],
                                             [-1, 0.5, 5]],
                                  cell=[10, 10, 10],
                                  pbc=True)

                    if wrap:
                        a.set_scaled_positions(a.get_scaled_positions() % 1.0)
                    a.set_calculator(LennardJones())
                    a.set_constraint(FixBondLength(0, 2))

                    d1 = a.get_distance(0, 2, mic=True)

                    FIRE(a, logfile=None).run(fmax=0.01)
                    e = a.get_potential_energy()
                    d2 = a.get_distance(0, 2, mic=True)
                    self.assertLess(abs(e - -2.034988), 1e-6)
                    self.assertLess(abs(d1 - d2), 1e-6)

    def test_mic(self):
        import ase
        import numpy as np

        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
                tol = 1e-9
                cell = np.array([[1., 0., 0.],
                                 [0.5, np.sqrt(3) / 2, 0.],
                                 [0., 0., 1.]]) * 10

                pos = np.dot(np.array([[0.0, 0.0, 0.0],
                                       [0.5, 0.5, 0.5],
                                       [0.2, 0.2, 0.2],
                                       [0.25, 0.5, 0.0]]), cell)

                a = Atoms('C4', pos, cell=cell, pbc=True)

                rpos = a.get_scaled_positions()

                # non-mic distance between atom 0 and 1
                d01F = np.linalg.norm(np.dot(rpos[1], cell))
                # mic distance between atom 0 (image [0,1,0]) and 1
                d01T = np.linalg.norm(np.dot(rpos[1] - np.array([0, 1, 0]), cell))
                d02F = np.linalg.norm(np.dot(rpos[2], cell))
                d02T = d02F
                # non-mic distance between atom 0 and 3
                d03F = np.linalg.norm(np.dot(rpos[3], cell))
                # mic distance between atom 0 (image [0,1,0]) and 3
                d03T = np.linalg.norm(np.dot(rpos[3] - np.array([0, 1, 0]), cell))

                # get_distance(mic=False)
                assert abs(a.get_distance(0, 1, mic=False) - d01F) < tol
                assert abs(a.get_distance(0, 2, mic=False) - d02F) < tol
                assert abs(a.get_distance(0, 3, mic=False) - d03F) < tol

                # get_distance(mic=True)
                assert abs(a.get_distance(0, 1, mic=True) - d01T) < tol
                assert abs(a.get_distance(0, 2, mic=True) - d02T) < tol
                assert abs(a.get_distance(0, 3, mic=True) - d03T) < tol

                # get_distance(mic=False, vector=True)
                assert all(abs(a.get_distance(0, 1, mic=False, vector=True)
                               - np.array([7.5, np.sqrt(18.75), 5.0])) < tol)
                assert all(abs(a.get_distance(0, 2, mic=False, vector=True)
                               - np.array([3., np.sqrt(3.), 2.0])) < tol)

                # get_distance(mic=True, vector=True)
                assert np.all(abs(a.get_distance(0, 1, mic=True, vector=True)
                                  - np.array([-2.5, np.sqrt(18.75), -5.0])) < tol)
                assert np.all(abs(a.get_distance(0, 2, mic=True, vector=True)
                                  - np.array([3., np.sqrt(3.), 2.0])) < tol)

                # get_all_distances(mic=False)
                all_dist = a.get_all_distances(mic=False)
                assert abs(all_dist[0, 1] - d01F) < tol
                assert abs(all_dist[0, 2] - d02F) < tol
                assert abs(all_dist[0, 3] - d03F) < tol
                assert all(abs(np.diagonal(all_dist)) < tol)

                # get_all_distances(mic=True)
                all_dist_mic = a.get_all_distances(mic=True)
                assert abs(all_dist_mic[0, 1] - d01T) < tol
                assert abs(all_dist_mic[0, 2] - d02T) < tol
                assert abs(all_dist_mic[0, 3] - d03T) < tol
                assert all(abs(np.diagonal(all_dist)) < tol)

                # get_distances(mic=False)
                for i in range(4):
                    assert all(abs(a.get_distances(i, [0, 1, 2, 3], mic=False) -
                                   all_dist[i]) < tol)

                # get_distances(mic=True)
                assert all(abs(a.get_distances(0, [0, 1, 2, 3], mic=True)
                               - all_dist_mic[0]) < tol)
                assert all(abs(a.get_distances(1, [0, 1, 2, 3], mic=True)
                               - all_dist_mic[1]) < tol)
                assert all(abs(a.get_distances(2, [0, 1, 2, 3], mic=True)
                               - all_dist_mic[2]) < tol)
                assert all(abs(a.get_distances(3, [0, 1, 2, 3], mic=True)
                               - all_dist_mic[3]) < tol)

                # get_distances(mic=False, vec=True)
                assert np.all(abs(a.get_distances(0, [0, 1, 2, 3], mic=False, vector=True)
                                  - np.array([a.get_distance(0, i, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)
                assert np.all(abs(a.get_distances(1, [0, 1, 2, 3], mic=False, vector=True)
                                  - np.array([a.get_distance(1, i, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)
                assert np.all(abs(a.get_distances(2, [0, 1, 2, 3], mic=False, vector=True)
                                  - np.array([a.get_distance(2, i, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)
                assert np.all(abs(a.get_distances(3, [0, 1, 2, 3], mic=False, vector=True)
                                  - np.array([a.get_distance(3, i, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)

                # get_distances(mic=True, vec=True)
                assert np.all(abs(a.get_distances(0, [0, 1, 2, 3], mic=True, vector=True)
                                  - np.array([a.get_distance(0, i, mic=True, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)
                assert np.all(abs(a.get_distances(1, [0, 1, 2, 3], mic=True, vector=True)
                                  - np.array([a.get_distance(1, i, mic=True, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)
                assert np.all(abs(a.get_distances(2, [0, 1, 2, 3], mic=True, vector=True)
                                  - np.array([a.get_distance(2, i, mic=True, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)
                assert np.all(abs(a.get_distances(3, [0, 1, 2, 3], mic=True, vector=True)
                                  - np.array([a.get_distance(3, i, mic=True, vector=True)
                                              for i in [0, 1, 2, 3]])) < tol)

                # set_distance
                a.set_distance(0, 1, 11., mic=False)
                assert abs(a.get_distance(0, 1, mic=False) - 11.) < tol
                assert abs(a.get_distance(0, 1, mic=True) - np.sqrt(46)) < tol

                # set_distance(mic=True)
                a.set_distance(0, 1, 3., mic=True)
                assert abs(a.get_distance(0, 1, mic=True) - 3.) < tol

    # def test_neighbor_kernel(self): # very slow (ignore for now)
    #     import ase.lattice.hexagonal
    #     from ase.build import bulk, molecule
    #     from ase.neighborlist import (mic, neighbor_list, primitive_neighbor_list,
    #                                   first_neighbors)
    #
    #     for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
    #         with self.subTest(msg = 'Atoms in ' + class_name):
    #             tol = 1e-7
    #             # two atoms
    #             a = Atoms('CC', positions=[[0.5, 0.5, 0.5], [1, 1, 1]], cell=[10, 10, 10],
    #                           pbc=True)
    #             i, j, d = neighbor_list("ijd", a, 1.1)
    #             assert (i == np.array([0, 1])).all()
    #             assert (j == np.array([1, 0])).all()
    #             assert np.abs(d - np.array([np.sqrt(3 / 4), np.sqrt(3 / 4)])).max() < tol
    #
    #             # test_neighbor_list
    #             for pbc in [True, False, [True, False, True]]:
    #                 a = Atoms('4001C', cell=[29, 29, 29])
    #                 a.set_scaled_positions(np.transpose([np.random.random(len(a)),
    #                                                      np.random.random(len(a)),
    #                                                      np.random.random(len(a))]))
    #                 j, dr, i, abs_dr, shift = neighbor_list("jDidS", a, 1.85)
    #
    #                 assert (np.bincount(i) == np.bincount(j)).all()
    #
    #                 r = a.get_positions()
    #                 dr_direct = mic(r[j] - r[i], a.cell)
    #                 assert np.abs(r[j] - r[i] + shift.dot(a.cell) - dr_direct).max() < tol
    #
    #                 abs_dr_from_dr = np.sqrt(np.sum(dr * dr, axis=1))
    #                 abs_dr_direct = np.sqrt(np.sum(dr_direct * dr_direct, axis=1))
    #
    #                 assert np.all(np.abs(abs_dr - abs_dr_from_dr) < 1e-12)
    #                 assert np.all(np.abs(abs_dr - abs_dr_direct) < 1e-12)
    #
    #                 assert np.all(np.abs(dr - dr_direct) < 1e-12)
    #
    #             # test_neighbor_list_atoms_outside_box
    #             for pbc in [True, False, [True, False, True]]:
    #                 a = Atoms('4001C', cell=[29, 29, 29])
    #                 a.set_scaled_positions(np.transpose([np.random.random(len(a)),
    #                                                      np.random.random(len(a)),
    #                                                      np.random.random(len(a))]))
    #                 a.set_pbc(pbc)
    #                 a.positions[100, :] += a.cell[0, :]
    #                 a.positions[200, :] += a.cell[1, :]
    #                 a.positions[300, :] += a.cell[2, :]
    #                 j, dr, i, abs_dr, shift = neighbor_list("jDidS", a, 1.85)
    #
    #                 assert (np.bincount(i) == np.bincount(j)).all()
    #
    #                 r = a.get_positions()
    #                 dr_direct = mic(r[j] - r[i], a.cell)
    #                 assert np.abs(r[j] - r[i] + shift.dot(a.cell) - dr_direct).max() < tol
    #
    #                 abs_dr_from_dr = np.sqrt(np.sum(dr * dr, axis=1))
    #                 abs_dr_direct = np.sqrt(np.sum(dr_direct * dr_direct, axis=1))
    #
    #                 assert np.all(np.abs(abs_dr - abs_dr_from_dr) < 1e-12)
    #                 assert np.all(np.abs(abs_dr - abs_dr_direct) < 1e-12)
    #
    #                 assert np.all(np.abs(dr - dr_direct) < 1e-12)
    #
    #             # test_small_cell
    #             a = Atoms('C', positions=[[0.5, 0.5, 0.5]], cell=[1, 1, 1],
    #                           pbc=True)
    #             i, j, dr, shift = neighbor_list("ijDS", a, 1.1)
    #             assert np.bincount(i)[0] == 6
    #             assert (dr == shift).all()
    #
    #             i, j = neighbor_list("ij", a, 1.5)
    #             assert np.bincount(i)[0] == 18
    #
    #             a.set_pbc(False)
    #             i = neighbor_list("i", a, 1.1)
    #             assert len(i) == 0
    #
    #             a.set_pbc([True, False, False])
    #             i = neighbor_list("i", a, 1.1)
    #             assert np.bincount(i)[0] == 2
    #
    #             a.set_pbc([True, False, True])
    #             i = neighbor_list("i", a, 1.1)
    #             assert np.bincount(i)[0] == 4
    #
    #             # test_out_of_cell_small_cell
    #             a = Atoms('CC', positions=[[0.5, 0.5, 0.5],
    #                                            [1.1, 0.5, 0.5]],
    #                           cell=[1, 1, 1], pbc=False)
    #             i1, j1, r1 = neighbor_list("ijd", a, 1.1)
    #             a.set_cell([2, 1, 1])
    #             i2, j2, r2 = neighbor_list("ijd", a, 1.1)
    #
    #             assert (i1 == i2).all()
    #             assert (j1 == j2).all()
    #             assert np.abs(r1 - r2).max() < tol
    #
    #             # test_out_of_cell_large_cell
    #             a = Atoms('CC', positions=[[9.5, 0.5, 0.5],
    #                                            [10.1, 0.5, 0.5]],
    #                           cell=[10, 10, 10], pbc=False)
    #             i1, j1, r1 = neighbor_list("ijd", a, 1.1)
    #             a.set_cell([20, 10, 10])
    #             i2, j2, r2 = neighbor_list("ijd", a, 1.1)
    #
    #             assert (i1 == i2).all()
    #             assert (j1 == j2).all()
    #             assert np.abs(r1 - r2).max() < tol
    #
    #             # test_hexagonal_cell
    #             for sx in range(3):
    #                 a = ase.lattice.hexagonal.Graphite('C', latticeconstant=(2.5, 10.0),
    #                                                    size=[sx + 1, sx + 1, 1])
    #                 i = neighbor_list("i", a, 1.85)
    #                 assert np.all(np.bincount(i) == 3)
    #
    #             # test_first_neighbors
    #             i = [1, 1, 1, 1, 3, 3, 3]
    #             assert (first_neighbors(5, i) == np.array([0, 0, 4, 4, 7, 7])).all()
    #             i = [0, 1, 2, 3, 4, 5]
    #             assert (first_neighbors(6, i) == np.array([0, 1, 2, 3, 4, 5, 6])).all()
    #
    #             # test_multiple_elements
    #             a = molecule('HCOOH')
    #             a.center(vacuum=5.0)
    #             i = neighbor_list("i", a, 1.85)
    #             assert (np.bincount(i) == np.array([2, 3, 1, 1, 1])).all()
    #
    #             cutoffs = {(1, 6): 1.2}
    #             i = neighbor_list("i", a, cutoffs)
    #             assert (np.bincount(i) == np.array([0, 1, 0, 0, 1])).all()
    #
    #             cutoffs = {(6, 8): 1.4}
    #             i = neighbor_list("i", a, cutoffs)
    #             assert (np.bincount(i) == np.array([1, 2, 1])).all()
    #
    #             cutoffs = {('H', 'C'): 1.2, (6, 8): 1.4}
    #             i = neighbor_list("i", a, cutoffs)
    #             assert (np.bincount(i) == np.array([1, 3, 1, 0, 1])).all()
    #
    #             cutoffs = [0.0, 0.9, 0.0, 0.5, 0.5]
    #             i = neighbor_list("i", a, cutoffs)
    #             assert (np.bincount(i) == np.array([0, 1, 0, 0, 1])).all()
    #
    #             cutoffs = [0.7, 0.9, 0.7, 0.5, 0.5]
    #             i = neighbor_list("i", a, cutoffs)
    #             assert (np.bincount(i) == np.array([2, 3, 1, 1, 1])).all()
    #
    #             # test_noncubic
    #             a = bulk("Al", cubic=False)
    #             i, j, d = neighbor_list("ijd", a, 3.1)
    #             assert (np.bincount(i) == np.array([12])).all()
    #             assert np.abs(d - [2.86378246] * 12).max() < tol
    #
    #             # test pbc
    #             nat = 10
    #             atoms = Atoms(numbers=range(nat),
    #                               cell=[(0.2, 1.2, 1.4),
    #                                     (1.4, 0.1, 1.6),
    #                                     (1.3, 2.0, -0.1)])
    #             atoms.set_scaled_positions(3 * np.random.random((nat, 3)) - 1)
    #
    #             for p1 in range(2):
    #                 for p2 in range(2):
    #                     for p3 in range(2):
    #                         atoms.set_pbc((p1, p2, p3))
    #                         i, j, d, D, S = neighbor_list("ijdDS", atoms, atoms.numbers * 0.2 + 0.5)
    #                         c = np.bincount(i, minlength=len(atoms))
    #                         atoms2 = atoms.repeat((p1 + 1, p2 + 1, p3 + 1))
    #                         i2, j2, d2, D2, S2 = neighbor_list("ijdDS", atoms2, atoms2.numbers * 0.2 + 0.5)
    #                         c2 = np.bincount(i2, minlength=len(atoms))
    #                         c2.shape = (-1, nat)
    #                         dd = d.sum() * (p1 + 1) * (p2 + 1) * (p3 + 1) - d2.sum()
    #                         dr = np.linalg.solve(atoms.cell.T, (atoms.positions[1] - atoms.positions[0]).T).T + np.array(
    #                             [0, 0, 3])
    #                         assert abs(dd) < 1e-10
    #                         assert not (c2 - c).any()
    #
    #             c = 0.0058
    #             i, j, d = primitive_neighbor_list('ijd',
    #                                               [True, True, True],
    #                                               np.eye(3) * 7.56,
    #                                               np.array([[0, 0, 0],
    #                                                         [0, 0, 0.99875]]),
    #                                               [c, c],
    #                                               self_interaction=False,
    #                                               use_scaled_positions=True)
    #             assert np.all(i == [0, 1])
    #             assert np.all(j == [1, 0])
    #             assert np.allclose(d, [0.00945, 0.00945])
    #
    #             # Empty atoms object
    #             i, D, d, j, S = neighbor_list("iDdjS", Atoms(), 1.0)
    #             assert i.dtype == np.int
    #             assert j.dtype == np.int
    #             assert d.dtype == np.float
    #             assert D.dtype == np.float
    #             assert S.dtype == np.int
    #             assert i.shape == (0,)
    #             assert j.shape == (0,)
    #             assert d.shape == (0,)
    #             assert D.shape == (0, 3)
    #             assert S.shape == (0, 3)
    #
    #             # Check that only a scalar (not a tuple) is returned if we request a single
    #             # argument.
    #             i = neighbor_list("i", Atoms(), 1.0)
    #             assert i.dtype == np.int
    #             assert i.shape == (0,)
    #
    def test_negativeindex(self):
        from ase.constraints import FixScaled
        for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
            with self.subTest(msg = 'Atoms in ' + class_name):
                a1 = Atoms(symbols='X2',
                           positions=[[0., 0., 0.], [2., 0., 0.], ],
                           cell=[[4., 0., 0.], [0., 4., 0.], [0., 0., 4.], ],
                           )

                fs1 = FixScaled(a1.get_cell(), -1, mask=(True, False, False))
                fs2 = FixScaled(a1.get_cell(), 1, mask=(False, True, False))

                a1.set_constraint([fs1, fs2])

                # reassigning using atoms.__getitem__
                a2 = a1[0:2]

                assert len(a1._constraints) == len(a2._constraints)

    # def cp2k_MD(self): # doesn't work with anything including Atoms!
    #     from ase import units
    #     from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    #     from ase.md.verlet import VelocityVerlet
    #     from ase.calculators.cp2k import CP2K
    #     for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
    #         with self.subTest(msg = 'Atoms in ' + class_name):
    #             calc = CP2K(label='test_H2_MD')
    #             positions = [(0, 0, 0), (0, 0, 0.7245595)]
    #             atoms = Atoms('HH', positions=positions, calculator=calc)
    #             atoms.center(vacuum=2.0)
    #
    #             # Run MD
    #             MaxwellBoltzmannDistribution(atoms, 0.5 * 300 * units.kB, force_temp=True)
    #             energy_start = atoms.get_potential_energy() + atoms.get_kinetic_energy()
    #             dyn = VelocityVerlet(atoms, 0.5 * units.fs)
    #             # def print_md():
    #             #    energy = atoms.get_potential_energy() + atoms.get_kinetic_energy()
    #             #    print("MD total-energy: %.10feV" %  energy)
    #             # dyn.attach(print_md, interval=1)
    #             dyn.run(20)
    #
    #             energy_end = atoms.get_potential_energy() + atoms.get_kinetic_energy()
    #
    #             assert energy_start - energy_end < 1e-4

    # def test_molecule(self):
    #     from ase.optimize import BFGS
    #     from ase.calculators.crystal import CRYSTAL
    #     for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
    #         with self.subTest(msg='Atoms in ' + class_name):
    #             with open('basis', 'w') as fd:
    #                 fd.write("""6 4
    #             0 0 6 2.0 1.0
    #              3048.0 0.001826
    #              456.4 0.01406
    #              103.7 0.06876
    #              29.23 0.2304
    #              9.349 0.4685
    #              3.189 0.3628
    #             0 1 2 4.0 1.0
    #              3.665 -0.3959 0.2365
    #              0.7705 1.216 0.8606
    #             0 1 1 0.0 1.0
    #              0.26 1.0 1.0
    #             0 3 1 0.0 1.0
    #              0.8 1.0
    #             """)
    #
    #             geom = Atoms('OHH',
    #                          positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])
    #
    #             geom.set_calculator(CRYSTAL(label='water',
    #                                         guess=True,
    #                                         basis='sto-3g',
    #                                         xc='PBE',
    #                                         otherkeys=['scfdir', 'anderson',
    #                                                    ['maxcycles', '500'],
    #                                                    ['toldee', '6'],
    #                                                    ['tolinteg', '7 7 7 7 14'],
    #                                                    ['fmixing', '90']]))
    #
    #             opt = BFGS(geom)
    #             opt.run(fmax=0.05)
    #
    #             final_energy = geom.get_potential_energy()
    #             assert abs(final_energy + 2047.34531091) < 1.0

    # def test_eon_multi_image_read(self):
    #     import tempfile
    #     import os
    #     import shutil
    #
    #     from numpy import array
    #     import ase.io
    #     for Atoms, class_name in zip(self.atoms_subclass_list, self.class_names):
    #         with self.subTest(msg='Atoms in ' + class_name):
    #             # Error tolerance.
    #             TOL = 1e-6
    #
    #             # A correct .con file.
    #             CON_FILE = """\
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                3.04618285858587523    2.15398224021542450    4.60622193000079427 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                3.36369427985916092    2.20887986058760699    4.61557342394151693 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                3.62116697589668135    2.40183843231018113    4.63674682635805180 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                3.83933949582109157    2.63825709178043821    4.65669727965894875 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                4.06162234392075128    2.87141228075409316    4.66819033729618926 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                4.28380999862949441    3.10482201783771883    4.65660558580467221 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                4.50188452903429326    3.34154720502221236    4.63664894132718874 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                4.75928919917819293    3.53496190773495211    4.61566200013953409 0   16
    #             0	Random Number Seed
    #             0	Time
    #             8.123222	5.744000	9.747867
    #             90.000000	90.000000	90.000000
    #             0 0
    #             0 0 0
    #             1
    #             17
    #             51.996100
    #             Cr
    #             Coordinates of Component 1
    #                1.01540277999999962    0.71799999999999997    1.01540277999999984 1    0
    #                3.04620834000000063    2.15399999999999991    1.01540277999999984 1    1
    #                3.04620834000000063    0.71799999999999997    3.04620834000000196 1    2
    #                1.01540277999999962    2.15399999999999991    3.04620834000000196 1    3
    #                1.01540277999999962    3.58999999999999986    1.01540277999999984 1    4
    #                3.04620834000000063    5.02599999999999980    1.01540277999999984 1    5
    #                3.04620834000000063    3.58999999999999986    3.04620834000000196 1    6
    #                1.01540277999999962    5.02599999999999980    3.04620834000000196 1    7
    #                5.07701389999999986    0.71799999999999997    1.01540277999999984 1    8
    #                7.10781945999998488    2.15399999999999991    1.01540277999999984 1    9
    #                7.10781945999998488    0.71799999999999997    3.04620834000000196 1   10
    #                5.07701389999999986    2.15399999999999991    3.04620834000000196 1   11
    #                5.07701389999999986    3.58999999999999986    1.01540277999999984 1   12
    #                7.10781945999998488    5.02599999999999980    1.01540277999999984 1   13
    #                7.10781945999998488    3.58999999999999986    3.04620834000000196 1   14
    #                5.07701389999999986    5.02599999999999980    3.04620834000000196 1   15
    #                5.07701160164306042    3.58998956883621734    4.60626159988447537 0   16
    #             """
    #
    #             # The corresponding data for the second to last image as an ASE Atoms object.
    #             data = ase.Atoms('Cr17', cell=array([[8.123222, 0, 0],
    #                                                  [0, 5.744000, 0],
    #                                                  [0, 0, 9.747867]]),
    #                              positions=array([[1.01540277999999962, 0.71799999999999997, 1.01540277999999984],
    #                                               [3.04620834000000063, 2.15399999999999991, 1.01540277999999984],
    #                                               [3.04620834000000063, 0.71799999999999997, 3.04620834000000196],
    #                                               [1.01540277999999962, 2.15399999999999991, 3.04620834000000196],
    #                                               [1.01540277999999962, 3.58999999999999986, 1.01540277999999984],
    #                                               [3.04620834000000063, 5.02599999999999980, 1.01540277999999984],
    #                                               [3.04620834000000063, 3.58999999999999986, 3.04620834000000196],
    #                                               [1.01540277999999962, 5.02599999999999980, 3.04620834000000196],
    #                                               [5.07701389999999986, 0.71799999999999997, 1.01540277999999984],
    #                                               [7.10781945999998488, 2.15399999999999991, 1.01540277999999984],
    #                                               [7.10781945999998488, 0.71799999999999997, 3.04620834000000196],
    #                                               [5.07701389999999986, 2.15399999999999991, 3.04620834000000196],
    #                                               [5.07701389999999986, 3.58999999999999986, 1.01540277999999984],
    #                                               [7.10781945999998488, 5.02599999999999980, 1.01540277999999984],
    #                                               [7.10781945999998488, 3.58999999999999986, 3.04620834000000196],
    #                                               [5.07701389999999986, 5.02599999999999980, 3.04620834000000196],
    #                                               [4.75928919917819293, 3.53496190773495211, 4.61566200013953409]]),
    #                              pbc=(True, True, True))
    #
    #             tempdir = tempfile.mkdtemp()
    #             try:
    #                 # First, write a correct .con file and try to read it.
    #                 con_file = os.path.join(tempdir, 'neb.con')
    #                 with open(con_file, 'w') as f:
    #                     f.write(CON_FILE)
    #                 images = ase.io.read(con_file, format='eon', index=':')
    #                 box = images[-2]
    #                 # Check cell vectors.
    #                 assert (abs(box.cell - data.cell)).sum() < TOL  # read: cell vector check
    #                 # Check atom positions.
    #                 # read: position check
    #                 assert (abs(box.positions - data.positions)).sum() < TOL
    #
    #                 # Now that we know that reading a .con file works, we will write
    #                 # one and read it back in.
    #                 out_file = os.path.join(tempdir, 'out.con')
    #                 ase.io.write(out_file, data, format='eon')
    #                 data2 = ase.io.read(out_file, format='eon')
    #                 # Check cell vectors.
    #                 # write: cell vector check
    #                 assert (abs(data2.cell - data.cell)).sum() < TOL
    #                 # Check atom positions.
    #                 # write: position check
    #                 assert (abs(data2.positions - data.positions)).sum() < TOL
    #             finally:
    #                 shutil.rmtree(tempdir)
    # def test_water(self):
    #     from ase.calculators.gaussian import Gaussian
    #     from ase.atoms import Atoms
    #     from ase.optimize.lbfgs import LBFGS
    #
    #     # First test to make sure Gaussian works
    #     calc = Gaussian(method='pbepbe', basis='sto-3g', force='force',
    #                     nproc=1, chk='water.chk', label='water')
    #     calc.clean()
    #
    #     water = Atoms('OHH',
    #                   positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)],
    #                   calculator=calc)
    #
    #     opt = LBFGS(water)
    #     opt.run(fmax=0.05)
    #
    #     forces = water.get_forces()
    #     energy = water.get_potential_energy()
    #     positions = water.get_positions()
    #
    #     # Then test the IO routines
    #     from ase.io import read
    #     water2 = read('water.log')
    #     forces2 = water2.get_forces()
    #     energy2 = water2.get_potential_energy()
    #     positions2 = water2.get_positions()
    #     # compare distances since positions are different in standard orientation
    #     dist = water.get_all_distances()
    #     dist2 = read('water.log', quantity='structures')[-1].get_all_distances()
    #
    #     assert abs(energy - energy2) < 1e-7
    #     assert abs(forces - forces2).max() < 1e-9
    #     assert abs(positions - positions2).max() < 1e-6
    #     assert abs(dist - dist2).max() < 1e-6
