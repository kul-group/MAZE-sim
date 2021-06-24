from maze.zeolite import PerfectZeolite, Zeolite, ClusterMaker
from maze.extra_framework_maker import ExtraFrameworkAnalyzer, ExtraFrameworkMaker
from ase import Atoms, db
from typing import Union, Tuple
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms, build
import numpy as np
from ase.visualize import view
import copy as copy
from ase.io import write, read
import itertools
from scipy.optimize import least_squares, minimize
import matplotlib.pyplot as plt
import pickle


class ForceFieldTrainer(object):

    def __init__(self, traj, initial_guess_param):
        self.traj = traj
        self.ref_forces = []
        self.bond = []
        self.angle = []
        self.good_atoms = []
        self.bad_atoms = []
        self.initial_guess_param = initial_guess_param
        self.ref_f_scalar = []

    @staticmethod
    def morse_force(r, a, e, rm):
        return 2 * a * e * np.exp(-a * (r - rm)) * (1 - np.exp(-a * (r - rm)))

    @staticmethod
    def harmonic(theta, k, theta_o):
        return 2 * k * (theta - theta_o)

    def get_properties(self):
        for atoms in self.traj:
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            EF_analyzer.get_extraframework()
            outlier = False
            """
            outlier = None
            for val in EF_analyzer.get_forces().values():
                if val[1] <= 0.5:
                    outlier = False
            """
            if outlier is False:
                self.bond.append(EF_analyzer.EF_bond_vectors)
                self.angle.append(EF_analyzer.EF_angles)
                self.ref_forces.append(EF_analyzer.get_forces())
                self.ref_f_scalar = self.get_all_vector_ref_forces()
                self.good_atoms.append(atoms)
            else:
                self.bad_atoms.append(atoms)   # todo: function handles bad atoms

    def get_angle_force_dir(self):
        ...

    def get_mul_predicted_forces(self, param):
        """
        f_harmonic = []
        for count, angle_dict in enumerate(self.angle):
            f_list = []
            for atom_tag, angle_list in angle_dict.items():
                f = []
                if 'O' in atom_tag:    # only considering Cu-O-Cu angle
                    for count, ang in enumerate(angle_list):
                        f += self.harmonic(ang, *param[3:5])


                f_list.append(f)

            f_harmonic.append(f_list)

        f_morse = []
        for count, bond_dict in enumerate(self.bond):
            f_list = []
            for atom_tag, dist_list in bond_dict.items():
                d_vec, d_mag, f = [], [], [0, 0, 0]
                for vec, mag in zip(dist_list[0], dist_list[1]):
                    d_vec.append(vec)
                    d_mag.append(mag)
                for count, vec in enumerate(d_vec):
                    f += self.morse_force(d_mag[count], *param[0:3]) * vec / d_mag[count]
                         # self.bonding_force(d_mag[count], *param[11:15]) * vec / d_mag[count]
                f_list.append(np.linalg.norm(f))
            assert len(f_list) == 3
            f_morse.append(f_list)
        """
        # use two sets of bond parameters
        f_morse = []
        for count, bond_dict in enumerate(self.bond):
            f_list = []
            for atom_tag, dist_list in bond_dict.items():
                d_vec, d_mag, f = [], [], [0, 0, 0]
                for vec, mag in zip(dist_list[0], dist_list[1]):
                    d_vec.append(vec)
                    d_mag.append(mag)
                if 'Cu' in atom_tag:
                    f += self.morse_force(d_mag[0], *param[3:6]) * d_vec[0] / d_mag[0]
                    f += self.morse_force(d_mag[1], *param[0:3]) * d_vec[1] / d_mag[1]
                if 'O' in atom_tag:
                    for count, vec in enumerate(d_vec):
                        f += self.morse_force(d_mag[count], *param[0:3]) * vec / d_mag[count]
                f_list.append(f)
            f_morse.append(f_list)

        f_all = np.array(f_morse)  #  + np.array(f_harmonic)
        return f_all

    def get_all_vector_ref_forces(self):
        vec_list = []
        for dicts in self.ref_forces:
            for atom_tag, force_list in dicts.items():
                vec_list.append(force_list[0])
        return vec_list

    def get_all_vector_predicted_forces(self, param):
        vec_list = []
        predicted_forces = self.get_mul_predicted_forces(param)
        for f_list in predicted_forces:
            for value in f_list:
                vec_list.append(value)
        return vec_list

    def get_residue(self, param):
        predicted_f_scalar = self.get_all_vector_predicted_forces(param)
        # print(np.array(predicted_f_scalar) - np.array(self.ref_f_scalar))
        return np.reshape(np.array(predicted_f_scalar) - np.array(self.ref_f_scalar), -1)[0::3]
        # residue = abs(np.array(predicted_f_scalar) - np.array(self.ref_f_scalar))
        # return sum(sum(residue))  # sum(residue[2::3])  # forces on O only

    def get_fitting_parameters(self):
        res = least_squares(self.get_residue, self.initial_guess_param, bounds=(0, np.inf))
            # minimize(self.get_residue, self.initial_guess_param, method='BFGS')
        #least_squares(self.get_residue, self.initial_guess_param)
        print(res)  # .success
        return res.x

    @property
    def optimized_params(self):
        return self.get_fitting_parameters()

    @property
    def optimized_forces(self):
        return self.get_all_vector_predicted_forces(self.optimized_params)

    def make_parity_plot(self):
        plt.figure()
        fig, ax = plt.subplots()
        x = [val[0] for val in self.ref_f_scalar]
        y = [val[0] for val in self.optimized_forces]
        plt.plot(x, y, 'o')
        plt.xlabel('DFT_force', fontsize=18)
        plt.ylabel('FF_force', fontsize=18)
        lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        # plt.title('fit: $\\epsilon:%4.2f, \\alpha:%4.2f, rm:%4.2f, k:%4.2f, \\theta_o:%4.2f, k:%4.2f, \\theta_o:%4.2f$' % tuple(self.optimized_params), fontsize=18)
        plt.show()

    def get_opt_dir(self, param, bonds, angles=None):
        # ignore angular terms for now
        f = [0, 0, 0]
        for count, vec in enumerate(bonds):
            f += self.morse_force(np.linalg.norm(vec), *param[0:3]) * vec / np.linalg.norm(vec)
        return f / np.linalg.norm(f)

    def optimize_geometry(self, atoms, params, coeff=0.005, cutoff_lim=0.1, cutoff=100):
        # eventually need three atoms, along three force directions, move with constant steps
        # now, moving the entire extra framework along the direction
        count = 1

        while abs(cutoff) >= cutoff_lim:
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            EF_analyzer.get_extraframework()
            EF_indices = EF_analyzer.EF_indices
            cutoff = EF_analyzer.get_predicted_forces(params)

            dir = self.get_opt_dir(params, EF_analyzer.EF_bond_vectors)
            original_pos = atoms.get_positions()[EF_indices]
            del atoms[EF_indices]

            atoms = atoms + Atoms('CuOCu', positions=original_pos + coeff * dir)
            print(abs(cutoff))
            count += 1
            if count == 300:
                break
        print('number of iterations =', count)


def test_trainer():
    traj = read('/Users/jiaweiguo/Box/ECH289_Project/MFI.traj', '0:2')
    ffTrainer = ForceFieldTrainer(traj, [1, 1, 2])
    ffTrainer.get_properties()
    # print(ffTrainer.get_residue([1, 1, 2, 1, 140]))

    print(ffTrainer.optimized_params)
    print(ffTrainer.optimized_forces)
    ffTrainer.make_parity_plot()


def test_optimizer():
    unopt_traj = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/MFI_CuOCu.traj', '0')
    params = [1, 1, 2, 1, 140]
    coeff = 0.05
    ffTrainer = ForceFieldTrainer(unopt_traj, params)
    # ffTrainer.optimize_geometry(unopt_traj, [1,1,1,1,1], coeff=0.005, cutoff_lim=0.1, cutoff=100)

    atoms = unopt_traj
    EF_analyzer = ExtraFrameworkAnalyzer(atoms)
    EF_analyzer.get_extraframework()
    EF_indices = EF_analyzer.EF_indices
    cutoff = EF_analyzer.get_predicted_forces(params)
    print(abs(cutoff))

    dir = ffTrainer.get_opt_dir(params, EF_analyzer.EF_bond_vectors)
    original_pos = atoms.get_positions()[EF_indices]
    del atoms[EF_indices]

    atoms = atoms + Atoms('OCuCu', positions=original_pos + coeff * dir)
    EF_analyzer = ExtraFrameworkAnalyzer(atoms)
    EF_analyzer.get_extraframework()
    cutoff = EF_analyzer.get_predicted_forces(params)
    print(abs(cutoff))


def test():
    traj = []
    F_list = []
    EF_index_list = []
    database = db.connect('/Users/jiaweiguo/Box/systematic_cuocu.db')
    for count, row in enumerate(database.select()):
        if count < 5:
            try:
                atoms = database.get_atoms(id=row.id)
                EF_analyzer = ExtraFrameworkAnalyzer(atoms)
                atoms = EF_analyzer.atoms
                index_list = EF_analyzer.get_extraframework_cluster()

                F = []
                EF_index_list.append(EF_analyzer.EF_indices)
                for index in EF_analyzer.EF_indices:
                    F.append(atoms.get_forces()[index])
                F_list.append(F)

                com_cell = np.matmul([0.5, 0.5, 0.5], atoms.get_cell())
                translate_vector = com_cell - atoms.get_positions()[EF_analyzer.centering_atom_index]
                new_atoms = atoms.translate(translate_vector)
                print(new_atoms.get_forces())
                atoms.wrap()

                assert len(index_list) == 13
                """
                # del atoms[[atom.index for atom in atoms if atom.index not in index_list]]
                # loss force information

                index_to_delete = [atom.index for atom in atoms if atom.index not in index_list]
                atoms = Zeolite(atoms)
                # print(atoms.get_forces())
                atoms = atoms.delete_atoms(index_to_delete)
                # print(atoms.get_forces())
                # print(atoms.get_forces())
                # new_atoms = atoms[[atom.index for atom in atoms if atom.index in index_list]]
                """

                traj.append(atoms)
            except:
                print('Error')

    print(EF_index_list)
    print(F_list)
    view(traj)

    my_data = {'CuOCu_index': EF_index_list, 'Forces': F_list}
    pickle.dump(my_data, open('CuOCu.p', 'wb'))
    write('/Users/jiaweiguo/Desktop/CuOCu_v2.traj ', traj)


if __name__ == '__main__':
    traj = read('/Users/jiaweiguo/Box/ECH289_Project/MFI.traj', '0::10')
    # read('/Users/jiaweiguo/Box/ECH289_Project/CHA.traj', ':')

    ini_param = [1, 1, 2, 1, 1, 4]
    ffTrainer = ForceFieldTrainer(traj, ini_param)
    ffTrainer.get_properties()
    # print(ffTrainer.bond)
    # print(ffTrainer.angle)
    # print(ffTrainer.ref_forces)
    # print(ffTrainer.get_mul_predicted_forces(ini_param))

    # print(ffTrainer.get_all_vector_ref_forces())
    # print(ffTrainer.get_all_scalar_predicted_forces(ini_param))

    # print(ffTrainer.get_residue(ini_param))

    # ffTrainer.get_fitting_parameters()

    # print(ffTrainer.optimized_params)
    # print(ffTrainer.optimized_forces)
    ffTrainer.make_parity_plot()

    # print(ffTrainer.get_residue(ffTrainer.optimized_params))


    # test_all()
    # test()
    # test_trainer()
    """
    traj = []
    F_list = []
    EF_index_list = []
    database = db.connect('/Users/jiaweiguo/Box/systematic_cuocu.db')
    for count, row in enumerate(database.select()):
        try:
            atoms = database.get_atoms(id=row.id)
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            atoms = EF_analyzer.atoms
            index_list = EF_analyzer.get_extraframework_cluster()

            F = []
            # EF_index_list.append(EF_analyzer.EF_indices)
            for index in EF_analyzer.EF_indices:
                F.append(atoms.get_forces()[index])
            F_list.append(F)
            
            com_cell = np.matmul([0.5, 0.5, 0.5], atoms.get_cell())
            translate_vector = com_cell - atoms.get_positions()[EF_analyzer.centering_atom_index]
            new_atoms = atoms.translate(translate_vector)
            atoms.wrap()
            
            assert len(index_list) == 13
            index_to_delete = [atom.index for atom in atoms if atom.index not in index_list]
            new_index = list((EF_analyzer.EF_indices - len(index_to_delete) * np.ones([1, 3]))[0])
            EF_index_list.append(new_index)
            assert new_index == [10, 11, 12]

            atoms = Zeolite(atoms)
            atoms = atoms.delete_atoms(index_to_delete)
            traj.append(atoms)
        except:
            print('Error')

    # print(EF_index_list)
    # print(F_list)

    my_data = {'Forces': F_list}
    pickle.dump(my_data, open('/Users/jiaweiguo/Desktop/CuOCu.p', 'wb'))
    write('/Users/jiaweiguo/Desktop/CuOCu_v2.traj', traj)

    traj = []
    F_list = []
    EF_index_list = []
    err_count = 0
    database = db.connect('/Users/jiaweiguo/Box/systematic_cuocu.db')
    for count, row in enumerate(database.select(id=1)):
        try:
            atoms = database.get_atoms(id=row.id)
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            atoms = EF_analyzer.atoms
            index_list = EF_analyzer.get_extraframework_cluster()

            assert len(index_list) == 13

            print(id(atoms.calc))
            atoms1 = copy.deepcopy(atoms)

            print(id(atoms1.calc))
            # F_list.append(atoms.calc.results['forces'][index_list])
            F_list.append(atoms.get_forces()[index_list])

            com_cell = np.matmul([0.5, 0.5, 0.5], atoms.get_cell())
            translate_vector = com_cell - atoms.get_positions()[EF_analyzer.centering_atom_index]
            new_atoms = atoms.translate(translate_vector)
            atoms.wrap()

            del atoms[[atom.index for atom in atoms if atom.index not in index_list]]
            traj.append(atoms)

            # print(atoms.calc.results['forces'])
        except:
            err_count += 1
            print('Error: ', err_count)
    # view(traj)

    my_data = {'Forces': F_list}
    pickle.dump(my_data, open('/Users/jiaweiguo/Desktop/CuOCu.p', 'wb'))
    write('/Users/jiaweiguo/Desktop/CuOCu_v2.traj', traj)
    """
