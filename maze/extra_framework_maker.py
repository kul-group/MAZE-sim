from maze.zeolite import PerfectZeolite, Zeolite
from ase import Atoms
from typing import Union, Tuple
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
import numpy as np
from ase.visualize import view
import copy as copy
from ase.io import write, read
import itertools


class ExtraFrameworkMaker(object):
    def __init__(self, iza_code=None):
        """ This is an extra-framework class
        :param iza_code: 3 letter code for the zeolite structures (IZA database)
        """
        if iza_code is not None:
            self.EFzeolite = Zeolite.make(iza_code)
        self.t_site_indices = {}
        self.traj_1Al = []
        self.traj_2Al = []
        self.count_all_Al_pairs = 0
        self.TM_list = ['Pt', 'Cu', 'Co', 'Pd', 'Fe', 'Cr', 'Rh', 'Ru']
        self.dict_t_sites_1Al_replaced = {}

    def make_extra_frameworks(self, print_statement=False):
        """
        :param print_statement:
        """
        self.replace_1Al()
        if print_statement is True:
            print('Single Al replacement is Done!')
        self.replace_2Al_unique_pairs()
        if print_statement is True:
            print('The second Al replacement is Done!')
        # add in Cu here

    def get_t_sites(self) -> dict:
        """ This function gets all unique T sites and corresponding atom indices for each T sites
        :return: dictionary mapping each unique T site tags with all atom indices with the same tag
        """
        t_site_indices = {}
        for site_name, value in self.EFzeolite.site_to_atom_indices.items():
            if 'T' in site_name:
                t_site_indices[site_name] = value
        return t_site_indices

    def replace_1Al(self):
        """ This function takes in a perfect zeolite and replace one Si atom with an Al for every T site
        :return: a dictionary of trajectories for each T site tags
        """
        self.t_site_indices = self.get_t_sites()
        for site_name, t_site in self.t_site_indices.items():
            traj_t_sites = []
            for index in t_site:
                pos = self.EFzeolite.get_positions()[index]
                new_zeo = self.EFzeolite.delete_atoms([index])
                new_zeo = new_zeo.add_atoms(Atoms('Al', positions=[pos]), 'Al')
                new_ztype = new_zeo.ztype + site_name + '->Al'
                new_zeo = Zeolite(new_zeo, ztype=new_ztype)
                self.traj_1Al.append(new_zeo)
                traj_t_sites.append(new_zeo)
            self.dict_t_sites_1Al_replaced[site_name] = traj_t_sites

    def replace_2Al_unique_pairs(self, cutoff_radius=9):
        """ This function makes the 2 Al replacement for all possible pairs (not limited to unique T-site pairs since
        even though the binding energies might be the same, the geometric properties, such as, Al-Al distance, are
        different). Replacement obeys the Lowenstein's rule, and Al pairs with distance greater than the "cutoff_radius"
        are ignored.
        :param cutoff_radius: replace the second Si site within some cutoff radius (9 Angstrom by default) around the
        first replacement site which is done using function "replace_1Al".
        :return:
        """
        done_indices = []

        for zeolite in self.traj_1Al:
            atoms = Zeolite(zeolite)  # already have neighbor_list
            ini_atoms = copy.copy(atoms)
            index_Al = [a.index for a in atoms if a.symbol == 'Al']

            if len(index_Al) != 1:
                print('Error: No more than 1 Al in the structure.')
                atoms[index_Al[len(index_Al) - 1]].symbol = 'Si'
                index_Al = [a.index for a in atoms if a.symbol == 'Al']

            neighboring_Si = []
            neigh_o_indices, offsets = atoms.neighbor_list.get_neighbors(index_Al[0])
            for each_neigh_o_index in neigh_o_indices:
                neigh_si_indices, offsets = atoms.neighbor_list.get_neighbors(each_neigh_o_index)
                for each_neigh_si_index in neigh_si_indices:
                    if each_neigh_si_index != index_Al[0]:
                        neighboring_Si.append(each_neigh_si_index)

            if len(neighboring_Si) != 4:
                continue
            else:
                Si_index = [a.index for a in atoms if a.symbol in ['Si']]
                for index in Si_index:
                    if index not in neighboring_Si and [index, index_Al] not in done_indices:
                        # first condition is the Lowenstein's rule
                        if 3.3 < atoms.get_distance(index_Al, index) < cutoff_radius:
                            T_site_name = ini_atoms.atom_indices_to_sites[index]
                            z_type_current = atoms.ztype
                            new_z_type = atoms.ztype + 'AND' + T_site_name + '->Al'
                            atoms = Zeolite(ini_atoms, ztype=new_z_type)
                            atoms[index].symbol = 'Al'
                            self.traj_2Al.append(atoms)
                            self.count_all_Al_pairs += 1
                            done_indices.append([index_Al, index])
                            done_indices.append([index, index_Al])

    @staticmethod
    def _get_close_framework_atom_positions(centering_index, atoms, radius_cutoff):
        """ This function detects all atoms on the input zeolite named "atoms" within certain "radius_cutoff" around the
        "centering_index" atom, and save all atoms positions into a list.
        :param centering_index: reference position of the extra-framework cluster to be inserted (currently placed with
        an S atom)
        :param atoms: the zeolite backbone
        :param radius_cutoff: 6 angstrom by default
        :return:
        """
        close_framework_atom_positions = []
        for count in range(len(atoms)):
            dist = atoms.get_distance(centering_index, count, mic=True)
            if dist <= radius_cutoff:
                close_framework_atom_positions.append(atoms[count].position)
        return close_framework_atom_positions

    def _get_reference_with_S_in_between_Al_pairs(self, atoms, cluster_radius):
        # TODO: write into smaller functions (check functions from the adsorbate class)
        """ This function finds the best position for the extra-framework cluster to be inserted (labeled as S for now).
        The reference S is placed in a way that it's not overlapping with the backbone atoms. The cluster size is taken
        into account from input "cluster_radius".
        :param atoms: zeolite framework with 2 Al atoms adjacent to each other
        :param cluster_radius: the approximated cluster radius
        :return:
        """
        Al_index = [a.index for a in atoms if a.symbol in ['Al']]
        mid_Al = 0.5 * (atoms.get_positions()[Al_index[0]] + atoms.get_positions()[Al_index[1]])
        atoms = atoms + Atoms('S', positions=[mid_Al])

        closest_distance = 1.5 + cluster_radius  # radius of Si atom = 1.5 Ang

        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        d_thres = 6.0  # might need to scale it with respect to the zeolite size
        close_framework_atom_positions = self._get_close_framework_atom_positions(index_S, atoms, d_thres)
        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]

        max_iter = 500
        num_iter = 0
        vec_Al = atoms.get_distance(Al_index[0], Al_index[1], mic=True, vector=True)
        u_Al = vec_Al / np.linalg.norm(vec_Al)  # unit direction along the Al-Al pair

        final_pos = []
        d_min = 0.0
        while num_iter <= max_iter:
            rand_vec = np.array([2.0 * np.random.random() - 1 for i in np.arange(3)])
            u_rand = rand_vec / np.linalg.norm(rand_vec)  # random unit vector
            u_dir = np.cross(u_Al, u_rand) / np.linalg.norm(np.cross(u_Al, u_rand))
            # some random vector perpendicular to the Al-Al vector
            u_dir = int(2 * np.round(np.random.random()) - 1) * u_dir
            # This takes a random number between 0 and 1, rounds it (0 or 1) and then shifts it to (-1 or 1)
            step_size = 0.5 * np.random.random_sample() * d_thres  # Pick step size randomly from 0 to half radius
            trial_pos = np.array(mid_Al + u_dir * step_size)

            distance_from_close_atoms = (close_framework_atom_positions - trial_pos)
            d_min = min(np.linalg.norm(distance_from_close_atoms, axis=1))

            num_iter += 1
            if closest_distance <= d_min:
                final_pos = trial_pos
                # closest_distance = d_min
                break

        cell_dimension = np.matmul(np.ones(3), atoms.get_cell())
        mic_constrained_pos = []
        for index, item in enumerate(final_pos):
            if item > cell_dimension[index]:
                item -= cell_dimension[index]
            elif item < 0:
                item += cell_dimension[index]
            mic_constrained_pos.append(item)

        atoms = atoms + Atoms('S', positions=[mic_constrained_pos])
        return atoms

    @staticmethod
    def _get_direction_of_insertion(atoms, index1, index2, index3):
        v1 = atoms.get_distance(index1, index2, vector=True, mic=True)
        v2 = atoms.get_distance(index1, index3, vector=True, mic=True)
        v = (v1 + v2) / np.linalg.norm(v1 + v2)
        return v

    def get_all_Z_TM(self, d_Z_TM, TM_type):
        """ # TODO: add description
        :param d_Z_TM: Z-TM distance with Z being the T sites on the zeolite framework and TM being extraframework atom
        to be inserted
        :return: a dictionary of structures for each T site name
        """
        dict_Z_TM = {}
        for site_name, all_zeo_with_same_T in self.dict_t_sites_1Al_replaced.items():
            atoms = copy.copy(all_zeo_with_same_T[0])
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            index_Al = [a.index for a in atoms if a.symbol == 'Al']
            indices, offsets = nl.get_neighbors(index_Al[0])
            assert len(indices) == 4

            traj = []
            pairs = list(itertools.combinations(indices, 2))
            for i, pair in enumerate(pairs):
                atoms = copy.copy(all_zeo_with_same_T[0])
                v = self._get_direction_of_insertion(atoms, index_Al[0], pair[0], pair[1])
                atoms = atoms + Atoms(TM_type, positions=[atoms[index_Al[0]].position] + v * d_Z_TM)
                traj.append(atoms)
            dict_Z_TM[site_name] = traj

        return dict_Z_TM

    def get_all_Bronsted_sites(self):
        traj_ZH = []
        done_oxygen = []
        for ref_atoms in self.traj_1Al:
            atoms = copy.copy(ref_atoms)
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            index_Al = [a.index for a in atoms if a.symbol == 'Al']
            indices, offsets = nl.get_neighbors(index_Al[0])
            assert len(indices) == 4

            for count, index in enumerate(indices):
                if index not in done_oxygen:
                    atoms = copy.copy(ref_atoms)
                    indices_neigh, offset_neigh = nl.get_neighbors(index)
                    assert len(indices_neigh) == 2  # oxygen index
                    i_neigh_T = [val for val in indices_neigh if val != index_Al][0]
                    assert atoms[i_neigh_T].symbol == 'Si'
                    assert atoms[index].symbol == 'O'
                    v = self._get_direction_of_insertion(atoms, index, i_neigh_T, index_Al)
                    coord_O = atoms.get_positions()[index]
                    traj_ZH.append(atoms.add_atoms(Atoms('H', positions=[coord_O + 0.97 * v]), 'H'))
                    done_oxygen.append(index)
        return traj_ZH

    @staticmethod
    def get_cluster_radius(EF_atoms):
        """ This function returns the averaged distance between atoms on the extra-framework cluster and the
        center-of-mass coordinate, which is used to represent/approximate the cluster size.
        :param EF_atoms: extra-framework cluster
        :return:
        """
        EF_center = EF_atoms.get_center_of_mass()
        all_position = EF_atoms.get_positions()
        radius_list = []
        for vector in abs(all_position - EF_center):
            radius_list.append(np.linalg.norm(vector))
        return np.mean(radius_list)

    def insert_ExtraFrameworkAtoms(self, my_zeolite, EF_atoms):
        """ This function takes in a zeolite backbone and an extra-framework cluster with the same cell dimensions as
        the zeolite. First, move the cluster center-of-mass to the reference position (indicated using an S atom). If
        there are more than one TMs in the cluster, the cluster is rotated so that the TM-TM vector is aligned with the
        Al-Al vector. Last, insert the cluster into "my_zeolite".
        :param my_zeolite: zeolite backbone with 2 Al atoms close to each other
        :param EF_atoms: the extra-framework cluster to be inserted in between the Al pair
        :return:
        """
        cluster_radius = self.get_cluster_radius(EF_atoms)
        atoms = self._get_reference_with_S_in_between_Al_pairs(my_zeolite, 0)  # cluster_radius = 0

        ## move S to the center of frame to aviod complicatio when inserting
        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        com_cell = np.matmul([0.5, 0.5, 0.5], atoms.get_cell())
        translate_vector = com_cell - atoms.get_positions()[index_S]
        atoms.translate(translate_vector)
        atoms.wrap()

        EF_atoms.set_cell(atoms.get_cell())
        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        ref_center = atoms.get_positions()[index_S]
        EF_atoms.translate(ref_center - EF_atoms.get_center_of_mass())
        EF_atoms.wrap()

        TM_index = [atom.index for atom in EF_atoms if atom.symbol in self.TM_list]
        Al_index = [atom.index for atom in atoms if atom.symbol == 'Al']

        if len(TM_index) == 2:
            ini_TM_vec = EF_atoms.get_positions()[TM_index[0]] - EF_atoms.get_positions()[TM_index[1]]
            final_vec = atoms.get_positions()[Al_index[0]] - atoms.get_positions()[Al_index[1]]
            EF_atoms.rotate(ini_TM_vec, final_vec, center=EF_atoms.get_center_of_mass())
            # ini_TM_vec is rotated into final_vec

        pos = EF_atoms.get_positions()
        atom_name = ''.join(EF_atoms.get_chemical_symbols())
        atoms = atoms + Atoms(atom_name, positions=pos, cell=atoms.get_cell())
        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]

        # Move the atoms back to the original position
        atoms.translate(-1 * translate_vector)
        atoms.wrap()
        return atoms
