from maze.zeolite import PerfectZeolite, Zeolite
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
from scipy.optimize import least_squares
import matplotlib.pyplot as plt


class ForceFieldTrainer(object):

    def __init__(self, traj, initial_guess_param):
        self.traj = traj
        self.ref_forces = []
        self.bond_list = []
        self.angle_list = []
        self.good_atoms = []
        self.bad_atoms = []
        self.initial_guess_param = initial_guess_param

    @staticmethod
    def morse_force(r, a, e, rm):
        return 2 * a * e * np.exp(-a * (r - rm)) * (1 - np.exp(-a * (r - rm)))

    @staticmethod
    def harmonic(theta, k, theta_o):
        return 2 * k * (theta - theta_o)

    def get_properties(self):
        for atoms in self.traj:
            try:
                EF_analyzer = ExtraFrameworkAnalyzer(atoms)
                EF_analyzer.get_extraframework_atoms()
                self.bond_list.append(EF_analyzer.EF_bond_vectors)
                self.angle_list.append(EF_analyzer.EF_angles)
                self.ref_forces.append(EF_analyzer.get_forces_on_centering_atom()[1])
                self.good_atoms.append(atoms)
            except:
                self.bad_atoms.append(atoms)   # todo: function handles bad atoms

    def get_mul_predicted_forces(self, param, bonds, angles):
        force_list = []
        for count, each_bond in enumerate(bonds):
            f = [0, 0, 0]
            for vec in each_bond:
                f += self.morse_force(np.linalg.norm(vec), *param[0:3]) * vec / np.linalg.norm(vec)
            force_list.append(np.linalg.norm(f) + self.harmonic(angles[count], *param[3:5]))
        return force_list

    def get_residue(self, param):
        return np.array(self.get_mul_predicted_forces(param, self.bond_list, self.angle_list)) - \
               np.array(self.ref_forces)

    def get_fitting_parameters(self):
        # TODO: FIND BETTER OPT ALGORITHM (DO IT FOR ECH289 PROJECT)
        # geo_properties = self.bond_list
        # geo_properties.extend((self.angle_list))
        res_lsq = least_squares(self.get_residue, self.initial_guess_param, bounds=(0, np.inf))
        return res_lsq.x

    @property
    def optimized_params(self):
        return self.get_fitting_parameters()

    @property
    def optimized_forces(self):
        return self.get_mul_predicted_forces(self.optimized_params, self.bond_list, self.angle_list)

    def make_parity_plot(self):
        plt.figure()
        fig, ax = plt.subplots()
        plt.plot(self.ref_forces, self.optimized_forces, 'o')
        plt.xlabel('DFT_force', fontsize=18)
        plt.ylabel('FF_force', fontsize=18)
        lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        plt.title('fit: $\\epsilon= %4.2f,  \\alpha= %4.2f,  rm= %4.2f, k= %4.2f,  \\theta_o= %4.2f$' %
                  tuple(self.optimized_params), fontsize=18)
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
            EF_analyzer.get_extraframework_atoms()
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
    traj = []
    database = db.connect('/Users/jiaweiguo/Box/systematic_cuocu.db')
    for row in database.select(min_energy_structure=True):
        traj.append(database.get_atoms(id=row.id))
    print(len(traj))

    ffTrainer = ForceFieldTrainer(traj, [1, 1, 2, 1, 140])
    ffTrainer.get_properties()
    # print(ffTrainer.get_residue([1, 1, 2, 1, 140]))

    print(ffTrainer.optimized_params)
    print(ffTrainer.optimized_forces)
    ffTrainer.make_parity_plot()


def get_unoptimized_zeo():
    sample_zeolite = 'MFI'
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))
    # write('/Users/jiaweiguo/Desktop/MFI_2Al.traj', EFzeolite.traj_2Al)

    EF_atoms = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
    MFI_2Al_traj = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/MFI_2Al.traj', ':')

    traj = []
    for zeolite in MFI_2Al_traj:
        EFzeolite = ExtraFrameworkMaker()
        traj.append(EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms))
    print('Extra-framework atoms insertion is done!')
    # write('/Users/jiaweiguo/Box/MAZE-sim-master/demos/MFI_CuOCu.traj', traj)


def test_optimizer():
    unopt_traj = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/MFI_CuOCu.traj', '0')
    params = [1, 1, 2, 1, 140]
    coeff = 0.05
    ffTrainer = ForceFieldTrainer(unopt_traj, params)
    # ffTrainer.optimize_geometry(unopt_traj, [1,1,1,1,1], coeff=0.005, cutoff_lim=0.1, cutoff=100)

    atoms = unopt_traj
    EF_analyzer = ExtraFrameworkAnalyzer(atoms)
    EF_analyzer.get_extraframework_atoms()
    EF_indices = EF_analyzer.EF_indices
    cutoff = EF_analyzer.get_predicted_forces(params)
    print(abs(cutoff))

    dir = ffTrainer.get_opt_dir(params, EF_analyzer.EF_bond_vectors)
    original_pos = atoms.get_positions()[EF_indices]
    del atoms[EF_indices]

    atoms = atoms + Atoms('OCuCu', positions=original_pos + coeff * dir)
    EF_analyzer = ExtraFrameworkAnalyzer(atoms)
    EF_analyzer.get_extraframework_atoms()
    cutoff = EF_analyzer.get_predicted_forces(params)
    print(abs(cutoff))




def test_all():
    """
    # get trained parameters
    traj = []
    database = db.connect('/Users/jiaweiguo/Box/systematic_cuocu.db')
    for row in database.select(min_energy_structure=True):
        traj.append(database.get_atoms(id=row.id))

    ffTrainer = ForceFieldTrainer(traj, [1, 1, 2, 1, 140])
    ffTrainer.get_properties()
    trained_params = ffTrainer.optimized_params
    """

    # pre-optimize
    unopt_traj = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/MFI_CuOCu.traj', '0:1')
    print(len(unopt_traj))
    dir_F = []
    for atoms in unopt_traj:
        EF_analyzer = ExtraFrameworkAnalyzer(atoms)
        EF_analyzer.get_extraframework_atoms()
        print(EF_analyzer.dict_EF_atoms)
        print(EF_analyzer.get_all_bonds())
        print(EF_analyzer.EF_indices)
        print(EF_analyzer.get_angle_at_centering_atom())
        # print(EF_analyzer.get_predicted_forces(trained_params))
        # dir_F.append(ffTrainer.get_opt_dir(trained_params, EF_analyzer.get_all_bonds()))
        """
        try:
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            EF_analyzer.get_extraframework_atoms()
            print(EF_analyzer.dict_EF_atoms)
            print(EF_analyzer.get_all_bonds())
            print(EF_analyzer.get_angle_at_centering_atom())
            print(EF_analyzer.get_forces_on_centering_atom())
            print(EF_analyzer.get_predicted_forces(trained_params))
            dir_F.append(ffTrainer.get_opt_dir(trained_params, EF_analyzer.get_all_bonds(),
                                  EF_analyzer.get_angle_at_centering_atom()))
        except:
            print('Error')
        """
    print(dir_F)
    # ADD Cu-Z terms


if __name__ == '__main__':
    # test_optimizer()
    atoms = read('/Users/jiaweiguo/Box/MAZE-sim-master/demos/MFI_CuOCu.traj', '4')
    EF_analyzer = ExtraFrameworkAnalyzer(atoms)
    EF_analyzer.get_extraframework_atoms()

