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
import time


class ForceFieldTrainer(object):

    def __init__(self, traj, initial_guess_param, scalar, tag):
        self.traj = traj
        self.ref_forces = []
        self.bond = []
        self.angle = []
        self.angle_dir = []
        self.good_atoms = []
        self.bad_atoms = []
        self.initial_guess_param = initial_guess_param
        self.ref_f_scalar = []
        self.optimized_forces = []
        self.optimized_params = []
        self.scalar = scalar
        self.tag = tag

    @staticmethod
    def morse_force(r, a, e, rm):
        return 2 * a * e * np.exp(-a * (r - rm)) * (1 - np.exp(-a * (r - rm)))

    @staticmethod
    def harmonic(theta, k, k1, theta_o):
        return 2 * k * (theta - theta_o) + 3 * k1 * (theta - theta_o) ** 2

    def get_properties(self):
        for atoms in self.traj:
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            EF_analyzer.get_extraframework()

            self.bond.append(EF_analyzer.EF_bond_vectors)
            self.angle.append(EF_analyzer.EF_angles)
            self.angle_dir.append(EF_analyzer.get_angle_force_dir())
            self.ref_forces.append(EF_analyzer.get_forces())
            self.good_atoms.append(atoms)
            # self.bad_atoms.append(atoms)   # todo: function handles bad atoms

        if self.scalar is True:
            self.ref_f_scalar = self.get_all_ref_forces()[1]
        else:
            self.ref_f_scalar = self.get_all_ref_forces()[0]

    def get_mul_predicted_forces(self, param):
        f_harmonic = []
        for count, angle_dict in enumerate(self.angle):
            f_list = []
            angle_dir = self.angle_dir[count]
            for index, (atom_tag, angle_list) in enumerate(angle_dict.items()):
                my_param = param[3:6]
                f = np.array(self.harmonic(angle_list[0], *my_param) * angle_dir[index]) / np.linalg.norm(angle_dir[index])
                f_list.append(f)
            if 'Cu' in self.tag:
                f_harmonic.append(f_list[0:2])  # only get forces on Cu
            if 'O' in self.tag:
                f_harmonic.append(f_list[2])  # only get forces on O

        f_morse = []
        for count, bond_dict in enumerate(self.bond):
            f_list = []
            for atom_tag, dist_list in bond_dict.items():
                d_vec, d_mag, f = [], [], [0, 0, 0]
                for vec, mag in zip(dist_list[0], dist_list[1]):
                    d_vec.append(vec)
                    d_mag.append(mag)
                if self.tag in atom_tag and 'Cu' in self.tag:
                    f += self.morse_force(d_mag[0], *param[0:3]) * d_vec[0] / d_mag[0]  # for Cu-O-Z bond
                    f += self.morse_force(d_mag[1], *param[6:9]) * d_vec[1] / d_mag[1]
                    # for Cu-O-Cu bond
                if self.tag in atom_tag and 'O' in self.tag:
                    for count, vec in enumerate(d_vec):
                        f += self.morse_force(d_mag[count], *param[0:3]) * vec / d_mag[count]
                f_list.append(f)
            if 'Cu' in self.tag:
                f_morse.append(f_list[0:2])  # only get forces on Cu
            if 'O' in self.tag:
                f_morse.append(f_list[2])  # only get forces on O

        f_all = np.array(f_morse) + np.array(f_harmonic)
        return f_all

    def get_all_ref_forces(self):
        vec_list = []
        sca_list = []
        for dicts in self.ref_forces:
            for atom_tag, force_list in dicts.items():
                if self.tag in atom_tag:
                    vec_list.append(force_list[0])
                    sca_list.append(force_list[1])
        return vec_list, sca_list

    def get_all_predicted_forces(self, param):
        vec_list = []
        sca_list = []
        predicted_forces = self.get_mul_predicted_forces(param)
        for f_list in predicted_forces:
            vec_list.append(f_list)
            sca_list.append(np.linalg.norm(f_list))
        return vec_list, sca_list

    def get_residue(self, param):
        if self.scalar is True:
            predicted_f = self.get_all_predicted_forces(param)[1]
        else:
            predicted_f = self.get_all_predicted_forces(param)[0]
        residue = np.reshape(np.array(np.reshape(predicted_f, [-1, 3])) - np.array(self.ref_f_scalar), -1)
        return np.mean(residue**2)

    def get_fitting_parameters(self):
        res = minimize(self.get_residue, self.initial_guess_param, method='nelder-mead')  # bounds=(0, np.inf)
            # minimize(self.get_residue, self.initial_guess_param, method='BFGS')
            #least_squares(self.get_residue, self.initial_guess_param)
        print(res.success)  # .success
        return res.x

    def get_optimized_properties(self):
        self.optimized_params = self.get_fitting_parameters()
        if self.scalar is True:
            self.optimized_forces = self.get_all_predicted_forces(self.optimized_params)[1]
        else:
            self.optimized_forces = self.get_all_predicted_forces(self.optimized_params)[0]

    def make_parity_plot(self):
        plt.figure()
        fig, ax = plt.subplots()
        self.get_optimized_properties()
        opt_f = np.array(np.reshape(self.optimized_forces, [-1, 3]))
        plt.plot(self.ref_f_scalar, opt_f, 'o')
        x = np.linspace(np.min(ax.get_xlim()), np.max(ax.get_xlim()), 100)
        coeff = np.polyfit(np.reshape(self.ref_f_scalar, -1), np.reshape(opt_f, -1), 1)
        y = np.polyval(coeff, x)
        plt.plot(x, y, 'r-')
        print('slope: ', coeff[0])
        plt.xlabel('DFT_force', fontsize=18)
        plt.ylabel('FF_force', fontsize=18)
        lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        plt.title('Force fitting on $%s$' %self.tag, fontsize=18)
        # plt.title('fit: $\\epsilon:%4.2f, \\alpha:%4.2f, rm:%4.2f, k:%4.2f, \\theta_o:%4.2f, k:%4.2f, \\theta_o:%4.2f$' % tuple(self.optimized_params), fontsize=18)
        plt.show()

    def get_MSE(self):
        rmspe = (np.sqrt(np.mean(np.square(np.reshape((np.array(self.ref_f_scalar) - np.array(self.optimized_forces)) / np.array(self.ref_f_scalar), -1))))) * 100
        return rmspe

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


def get_good_initial_traj(traj, tag, round_num, f_lim):
    good_traj = []
    x_list = []
    y_list = []
    z_list = []
    for atoms in traj:
        EF_analyzer = ExtraFrameworkAnalyzer(atoms)
        EF_analyzer.get_extraframework()
        val = [f[0] for f in EF_analyzer.get_forces().values()]
        if tag == 'Cu':
            val = val[0]
        if tag == 'O':
            val = val[2]
        if np.round(val[0], round_num) not in x_list and np.round(val[1], round_num) not in y_list and \
                np.round(val[2], round_num) not in z_list and abs(val[0]) < f_lim \
                and abs(val[1]) < f_lim and abs(val[2]) < f_lim:
            x_list.append(np.round(val[0], round_num))
            y_list.append(np.round(val[1], round_num))
            z_list.append(np.round(val[2], round_num))
            good_traj.append(atoms)
    return good_traj


if __name__ == '__main__':
    """
    traj = read('/Users/jiaweiguo/Box/MFI_all.traj', ':')
    # read('/Users/jiaweiguo/Box/MFI_all_Cu.traj', ':')
    # read('/Users/jiaweiguo/Box/ECH289_Project/MFI.traj', ':')
    # traj = read('/Users/jiaweiguo/Box/MFI_minE_Cu.traj', ':')
    # read('/Users/jiaweiguo/Box/MFI_all.traj', ':')
    # read('/Users/jiaweiguo/Box/ECH289_Project/MFI.traj', ':')
    # read('/Users/jiaweiguo/Box/ECH289_Project/CHA.traj', ':')
    tic = time.perf_counter()
    good_traj = get_good_initial_traj(traj, tag='O', round_num=1, f_lim=1.5)
    toc = time.perf_counter()
    print(f"Function terminated in {toc - tic:0.4f} seconds")
    write('/Users/jiaweiguo/Box/MFI_all_O.traj', good_traj)
    print(len(good_traj))
    """
    traj = read('/Users/jiaweiguo/Box/MFI_all_O.traj', ':')

    my_tag = 'O'

    if my_tag is 'O':
        ini_param = [1, 1, 2, 1, 1, 4]
    else:
        ini_param = [1, 1, 2, 1, 1, 4, 1, 1, 2]
    ffTrainer = ForceFieldTrainer(traj, ini_param, scalar=False, tag=my_tag)
    ffTrainer.get_properties()
    # print(ffTrainer.bond)
    # print(ffTrainer.angle)
    # print(ffTrainer.ref_forces)
    # print(ffTrainer.get_mul_predicted_forces(ini_param))

    # print(ffTrainer.get_all_ref_forces()[1])
    # print(ffTrainer.get_all_predicted_forces(ini_param))

    # print(ffTrainer.get_residue(ini_param))

    # ffTrainer.get_fitting_parameters()

    ffTrainer.make_parity_plot()
    print(ffTrainer.optimized_params)
    # print(ffTrainer.get_MSE())
    # print(ffTrainer.optimized_forces)

    # print(ffTrainer.get_residue(ffTrainer.optimized_params))
    # """
