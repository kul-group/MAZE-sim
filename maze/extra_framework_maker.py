from maze.zeolite import PerfectZeolite, Zeolite
from maze.io_zeolite import read_vasp
from ase import Atoms, db
from typing import Union, Tuple
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
import numpy as np
from ase.visualize import view
import copy as copy
from ase.io import write, read
import itertools
from scipy.optimize import least_squares


class ExtraFrameworkMaker(object):

    def __init__(self, iza_code=None, optimized_zeolite_path=None, user_input_path=None):
        """ This is an extra-framework class
        :param iza_code: 3 letter code for the zeolite structures (IZA database)
        """
        if iza_code is not None:
            self.EFzeolite = Zeolite.make(iza_code)
        if optimized_zeolite_path is not None:
            read_vasp(optimized_zeolite_path, Zeolite.make(iza_code))
        if user_input_path is not None:
            # self.EFzeolite = read(user_input_path, '0')
            self.EFzeolite = Zeolite(PerfectZeolite.build_from_cif_with_labels(filepath=user_input_path))

        self.t_site_indices = {}
        self.t_site_indices_count = {}
        self.traj_1Al = []
        self.traj_2Al = []
        self.count_all_Al_pairs = 0
        self.TM_list = ['Pt', 'Cu', 'Co', 'Pd', 'Fe', 'Cr', 'Rh', 'Ru']
        self.dict_1Al_replaced = {}
        self.dict_2Al_replaced = {}
        self.T_site_pair = []
        self.T_index_pair = []

    def make_extra_frameworks(self, replace_1Al=False, replace_2Al=False, print_statement=False):
        """
        :param replace_1Al:
        :param replace_2Al:
        :param print_statement:
        """
        if replace_1Al is True:
            self.replace_1Al()
            if print_statement is True:
                print('Single Al replacement is Done!')

        if replace_2Al is True:
            self.replace_2Al_unique_pairs()
            for count_num, item in enumerate(self.T_site_pair):
                key_tag = '_'.join(item) + '_' + '_'.join(map(str, self.T_index_pair[count_num]))
                self.dict_2Al_replaced[key_tag] = self.traj_2Al[count_num]
            if print_statement is True:
                print('The second Al replacement is Done!')
        # add in Cu here

    def get_t_sites(self):
        """ This function gets all unique T sites and corresponding atom indices for each T sites
        :return: dictionary mapping each unique T site tags with all atom indices with the same tag
        """
        for site_name, value in self.EFzeolite.site_to_atom_indices.items():
            if 'T' in site_name:
                self.t_site_indices[site_name] = value
                self.t_site_indices_count[site_name] = len(value)

    def replace_1Al(self):
        """ This function takes in a perfect zeolite and replace one Si atom with an Al for every T site
        :return: a dictionary of trajectories for each T site tags
        """
        self.get_t_sites()
        for site_name, t_site in self.t_site_indices.items():
            traj_t_sites = []
            for index in t_site:
                pos = self.EFzeolite.get_positions()[index]
                new_zeo = self.EFzeolite.delete_atoms([index])
                new_zeo = new_zeo.add_atoms(Atoms('Al', positions=[pos]), 'Al')
                # new_zeo = self.EFzeolite
                # new_zeo[index].symbol = 'Al'
                new_ztype = new_zeo.ztype + site_name + '->Al'
                new_zeo = Zeolite(new_zeo, ztype=new_ztype)
                self.traj_1Al.append(new_zeo)
                traj_t_sites.append(new_zeo)
            self.dict_1Al_replaced[site_name] = traj_t_sites

    def get_T_site_from_atom_index(self, index):
        return [k for k, v in self.t_site_indices.items() if index in v]

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
                            self.T_site_pair.append([self.get_T_site_from_atom_index(index_Al[0])[0], T_site_name])
                            self.T_index_pair.append([index_Al[0], index])
                            z_type_current = atoms.ztype
                            new_z_type = atoms.ztype + 'AND' + T_site_name + '->Al'
                            atoms = Zeolite(ini_atoms, ztype=new_z_type)
                            atoms[index].symbol = 'Al'
                            self.traj_2Al.append(atoms)
                            self.count_all_Al_pairs += 1
                            done_indices.append([index_Al, index])
                            done_indices.append([index, index_Al])

    @staticmethod
    def _get_direction_of_insertion(atoms, index1, index2, index3):
        v1 = atoms.get_distance(index1, index2, vector=True, mic=True)
        v2 = atoms.get_distance(index1, index3, vector=True, mic=True)
        v = (v1 + v2) / np.linalg.norm(v1 + v2)
        return v

    def get_all_Z_TM(self, d_Z_TM, TM_type):
        """
        :param d_Z_TM: Z-TM distance with Z being the T sites on the zeolite framework and TM being extraframework atom
        to be inserted
        :return: a dictionary of structures for each T site name
        """
        dict_Z_TM = {}
        for site_name, all_zeo_with_same_T in self.dict_1Al_replaced.items():
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

    def get_Z_TM(self, atoms, d_Z_TM, TM_type):

        nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
        nl.update(atoms)
        index_Al = [a.index for a in atoms if a.symbol == 'Al']
        indices, offsets = nl.get_neighbors(index_Al[0])
        assert len(indices) == 4

        dict_Z_TM = {}
        original_atoms = copy.copy(atoms)
        pairs = list(itertools.combinations(indices, 2))
        for i, pair in enumerate(pairs):
            atoms = copy.copy(original_atoms)
            v = self._get_direction_of_insertion(atoms, index_Al[0], pair[0], pair[1])
            atoms = atoms + Atoms(TM_type, positions=[atoms[index_Al[0]].position] + v * d_Z_TM)
            atoms.wrap()
            key_tag = 'O' + str(pair[0]) + '_O' + str(pair[1])
            dict_Z_TM[key_tag] = atoms

        return dict_Z_TM

    def _insert_H(self, atoms, index):
        nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
        nl.update(atoms)
        index_Al = [a.index for a in atoms if a.symbol == 'Al']
        indices_neigh, offset_neigh = nl.get_neighbors(index)
        assert len(indices_neigh) == 2
        i_neigh_T = [val for val in indices_neigh if val not in index_Al][0]
        assert atoms[i_neigh_T].symbol == 'Si'
        assert atoms[index].symbol == 'O'
        v = self._get_direction_of_insertion(atoms, index, i_neigh_T, index_Al)
        coord_O = atoms.get_positions()[index]
        new_atoms = atoms + Atoms('H', positions=[coord_O + 0.97 * v])
        return new_atoms

    def get_all_Bronsted_sites(self, case_1Al=False, case_2Al=False, my_dict=None):
        """
        This function samples all Bronsted sites for each T sites names (4 Bronsted sites for each T site).
        No double counting within each sites_name ('T1', 'T2', 'T3', etc), but do have overlaps among different
        site_names.
        """

        if case_1Al is True:
            my_dict = copy.copy(self.dict_1Al_replaced)
        if case_2Al is True:
            my_dict = copy.copy(self.dict_2Al_replaced)

        dict_ZH = {}
        for site_name, all_zeo_with_same_T in my_dict.items():
            if case_2Al is True:
                atoms = copy.copy(all_zeo_with_same_T)
            else:
                atoms = copy.copy(all_zeo_with_same_T[0])

            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            index_Al = [a.index for a in atoms if a.symbol == 'Al']
            indices = []
            for each_Al in index_Al:
                indices.append(list(nl.get_neighbors(each_Al)[0]))

            traj = []
            if case_2Al is True:
                assert len(indices) == 2
                all_pairs = []
                for count, value1 in enumerate(indices[0]):
                    all_pairs.extend([value1, value2] for value2 in indices[1])
                assert len(all_pairs) == 16
                for index_pair in all_pairs:
                    new_atoms = copy.copy(atoms)
                    for index in index_pair:
                        new_atoms = self._insert_H(new_atoms, index)
                    traj.append(new_atoms)
            else:
                for count, index in enumerate(indices[0]):
                    new_atoms = self._insert_H(atoms, index)
                    traj.append(new_atoms)

            dict_ZH[site_name] = traj
        return dict_ZH

    def get_Bronsted_sites(self, atoms):
        # for each input structure
        nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
        nl.update(atoms)
        index_Al = [a.index for a in atoms if a.symbol == 'Al']
        indices = []
        for each_Al in index_Al:
            indices.append(list(nl.get_neighbors(each_Al)[0]))

        dict_H = {}
        if len(index_Al) == 2:
            all_pairs = []
            for count, value1 in enumerate(indices[0]):
                all_pairs.extend([value1, value2] for value2 in indices[1])
            assert len(all_pairs) == 16
            for index_pair in all_pairs:
                new_atoms = copy.copy(atoms)
                for index in index_pair:
                    new_atoms = self._insert_H(new_atoms, index)
                key_tag = 'O' + str(index_pair[0]) + '_O' + str(index_pair[1])
                dict_H[key_tag] = new_atoms
        else:
            for count, index in enumerate(indices[0]):
                new_atoms = self._insert_H(atoms, index)
                key_tag = 'O' + str(index)
                dict_H[key_tag] = new_atoms
        return dict_H

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

    def _get_random_dir(self, atoms):
        Al_index = [a.index for a in atoms if a.symbol in ['Al']]
        vec_Al = atoms.get_distance(Al_index[0], Al_index[1], mic=True, vector=True)
        u_Al = vec_Al / np.linalg.norm(vec_Al)  # unit direction along the Al-Al pair

        rand_vec = np.array([2.0 * np.random.random() - 1 for i in np.arange(3)])
        u_rand = rand_vec / np.linalg.norm(rand_vec)  # random unit vector
        u_dir = np.cross(u_Al, u_rand) / np.linalg.norm(np.cross(u_Al, u_rand))
        # some random vector perpendicular to the Al-Al vector
        u_dir = int(2 * np.round(np.random.random()) - 1) * u_dir
        # This takes a random number between 0 and 1, rounds it (0 or 1) and then shifts it to (-1 or 1)
        return u_dir

    def _get_reference_with_S_in_between_Al_pairs(self, atoms, cluster_radius):
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

        index_S = [a.index for a in atoms if a.symbol == 'S'][0]
        d_thres = 6.0  # might need to scale it with respect to the zeolite size
        close_framework_atom_positions = self._get_close_framework_atom_positions(index_S, atoms, d_thres)
        del atoms[[atom.index for atom in atoms if atom.symbol == 'S']]

        max_iter, num_iter, d_min, closest_distance = 1000, 0, 0, 1.5 + cluster_radius  # radius of Si atom = 1.5 Ang
        final_pos = []
        while num_iter <= max_iter:
            u_dir = self._get_random_dir(atoms)
            step_size = d_thres * np.random.random_sample()  # Pick step size randomly from 0 to closest_distance
            trial_pos = np.array(mid_Al + u_dir * step_size)

            distance_from_close_atoms = (close_framework_atom_positions - trial_pos)
            d_min = min(np.linalg.norm(distance_from_close_atoms, axis=1))
            print(d_min)

            num_iter += 1
            if closest_distance <= d_min:
                final_pos = trial_pos
                # closest_distance = d_min
                break
        print(final_pos)

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


class ExtraFrameworkAnalyzer(object):

    def __init__(self, atoms):
        self.TM_list = ['Pt', 'Cu', 'Co', 'Pd', 'Fe', 'Cr', 'Rh', 'Ru']
        self.dict_EF = {}
        self.atoms = PerfectZeolite(atoms)

    @property
    def EF_bond_vectors(self):
        return self.get_all_bonds()

    @property
    def EF_angles(self):
        return self.get_all_angles()

    """
    def get_extraframework_atoms(self):
        #TODO: MORE GENERAL EXTRA FRAMEWORK DETECTION (CUFRRENTLY LIMITED TO TM-O-TM)

        index_EF_TM = [a.index for a in self.atoms if a.symbol in self.TM_list]
        index_Al = [a.index for a in self.atoms if a.symbol == 'Al']
        assert len(index_EF_TM) == 2
        assert len(index_Al) == 2

        # self.atoms.update_nl(1.2) need larger cutoff for tracking ZOCu oxygens

        TM_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_EF_TM[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_EF_TM[1])[0]))
        Al_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_Al[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_Al[1])[0]))

        # print(TM_neigh_list, Al_neigh_list)

        # This is wrong! Not always return desired O index!
        centering_o = [[x for x in TM_neigh_list if list(TM_neigh_list).count(x) > 1][0]]
        # print(centering_o)

        o_between_T_Cu = [val for val in TM_neigh_list if val in Al_neigh_list and self.atoms[val].symbol == 'O']
        # print(o_between_T_Cu)

        self.centering_atom_index = centering_o[0]
        assert len(centering_o) == 1
        assert self.atoms[centering_o].symbols == 'O'

        EF_indices = [index_EF_TM]
        EF_indices.extend(centering_o)
        EF_symbols = [self.atoms[index_EF_TM[0]].symbol]
        EF_symbols.extend('O')

        self.EF_indices = list(centering_o)
        self.EF_indices.extend([value for value in index_EF_TM])

        for count, index in enumerate(EF_indices):
            self.dict_EF_atoms[EF_symbols[count]] = index

        self.o_between_T_Cu = o_between_T_Cu
        # self.dict_EF_atoms['OZ'] = self.o_between_T_Cu
        
    def get_extraframework_cluster(self):
        # extraframework atoms, 2 Al and surrounding 8 O
        index_EF_TM = [a.index for a in self.atoms if a.symbol in self.TM_list]
        index_Al = [a.index for a in self.atoms if a.symbol == 'Al']
        assert len(index_EF_TM) == 2
        assert len(index_Al) == 2

        TM_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_EF_TM[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_EF_TM[1])[0]))
        Al_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_Al[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_Al[1])[0]))
        centering_o = [[x for x in TM_neigh_list if list(TM_neigh_list).count(x) > 1][0]]
        self.centering_atom_index = centering_o[0]
        assert len(centering_o) == 1
        assert self.atoms[centering_o].symbols == 'O'

        self.EF_indices = list(centering_o)
        self.EF_indices.extend([value for value in index_EF_TM])

        return np.unique(list(Al_neigh_list) + centering_o + index_Al + index_EF_TM)
    """

    def get_O_index_between_atoms(self, index_1, index_2, scale=3.0, O_count=2):
        # find the closest O in between two atoms since nl of ASE is so annoying

        self.atoms.update_nl(scale)
        nl_1 = self.atoms.neighbor_list.get_neighbors(index_1)[0]
        nl_2 = self.atoms.neighbor_list.get_neighbors(index_2)[0]
        index_O = [val for val in nl_1 if val in nl_2 and self.atoms[val].symbol == 'O']

        index_list = []
        dist_list = []
        for index in index_O:
            index_list.append(index)
            dist_list.append(0.5 * (self.atoms.get_distance(index_1, index, mic=True) + self.atoms.get_distance(index_2, index, mic=True)))

        unsorted_dist_list = copy.copy(dist_list)
        dist_list.sort()
        closed_O_index = []
        if O_count == 1:
            for index, element in enumerate(unsorted_dist_list):
                if element == dist_list[0]:
                    closed_O_index.append(index_list[index])
        else:
            for index, element in enumerate(unsorted_dist_list):
                if element == dist_list[0]:
                    closed_O_index.append(index_list[index])
                if element == dist_list[1]:
                    closed_O_index.append(index_list[index])
        return closed_O_index

    def get_extraframework(self):
        index_Al = [a.index for a in self.atoms if a.symbol == 'Al']
        index_Cu = [a.index for a in self.atoms if a.symbol == 'Cu']
        index_Al1, index_Al2 = index_Al[0], index_Al[1]
        if self.atoms.get_distance(index_Al1, index_Cu[0], mic=True) < self.atoms.get_distance(index_Al1, index_Cu[1], mic=True):
            index_Cu1, index_Cu2 = index_Cu[0], index_Cu[1]
        else:
            index_Cu1, index_Cu2 = index_Cu[1], index_Cu[0]
        centering_O = [288]  # [108]  # self.get_O_index_between_atoms(index_Cu1, index_Cu2, scale=0.85, O_count=1)
        Cu1_O_neigh = self.get_O_index_between_atoms(index_Al1, index_Cu1)
        Cu2_O_neigh = self.get_O_index_between_atoms(index_Al2, index_Cu2)
        self.dict_EF['Cu1'] = [index_Cu1, [index_Al1]+centering_O]  # [index_Cu1, Cu1_O_neigh+centering_O]
        self.dict_EF['Cu2'] = [index_Cu2, [index_Al2]+centering_O]  # [index_Cu2, Cu2_O_neigh+centering_O]
        self.dict_EF['O'] = [centering_O[0], [index_Cu1, index_Cu2]]

    def get_all_bonds(self):
        dict_EF_bonds = {}
        for atom_tag, index_list in self.dict_EF.items():
            atom_index = index_list[0]
            neighbor_index = index_list[1]
            d_vec = []
            d_mag = []
            for index in neighbor_index:
                d_vec.append(self.atoms.get_distance(index, atom_index, mic=True, vector=True))
                d_mag.append(self.atoms.get_distance(index, atom_index, mic=True, vector=False))
            dict_EF_bonds[atom_tag] = [d_vec, d_mag]
        return dict_EF_bonds

    def get_all_angles(self):
        dict_EF_angles = {}
        for atom_tag, index_list in self.dict_EF.items():
            atom_index = index_list[0]
            """
            neighbor_index = index_list[1][0:2]
            if 'Cu' in atom_tag:
                angle = []
                for index in neighbor_index:
                    angle.append(self.atoms.get_angle(index, atom_index, index_list[1][2], mic=True) / 180 * np.pi)  # O, Cu, O
            else:
                angle = [self.atoms.get_angle(neighbor_index[0], atom_index, neighbor_index[1], mic=True) / 180 * np.pi]
            """
            neighbor_index = index_list[1]
            angle = [self.atoms.get_angle(neighbor_index[0], atom_index, neighbor_index[1], mic=True) / 180 * np.pi]
            dict_EF_angles[atom_tag] = angle
        return dict_EF_angles

    def get_angle_force_dir(self):
        angle_dir = []
        dict_EF_bonds = self.get_all_bonds()
        for atom_tag, val in dict_EF_bonds.items():
            vector = val[0]
            angle_dir.append(-0.5 * (vector[0] + vector[1]))
        return angle_dir

    def get_forces(self):
        dict_EF_forces = {}
        for atom_tag, index_list in self.dict_EF.items():
            atom_index = index_list[0]
            f_vec = self.atoms.calc.results['forces'][atom_index]  # self.atoms.get_forces()[atom_index]
            f_mag = np.linalg.norm(f_vec)
            dict_EF_forces[atom_tag] = [f_vec, f_mag]
        return dict_EF_forces

    """
    @staticmethod
    def bonding_force(r, K1, K2, K3, R):
        return 2 * K1 * (r - R) + 3 * K2 * (r - R) ** 2 + 4 * K3 * (r - R) ** 3

    @staticmethod
    def morse_force(r, a, e, rm):
        return 2 * a * e * np.exp(-a * (r - rm)) * (1 - np.exp(-a * (r - rm)))

    @staticmethod
    def harmonic(theta, k, theta_o):
        return 2 * k * (theta - theta_o)

    def get_predicted_forces(self, param):
        # param = [a, e, rm, k_Cu, theta_o_Cu, k_O, theta_o_O, K, R]

        dict_EF_predicted_forces = {}
        dict_EF_bonds = self.get_all_bonds()
        dict_EF_angles = self.get_all_angles()

        dict_morse = {}
        for atom_tag, dist_list in dict_EF_bonds.items():
            d_vec, d_mag, f = [], [], [0, 0, 0]
            for vec, mag in zip(dist_list[0], dist_list[1]):
                d_vec.append(vec)
                d_mag.append(mag)
            for count, vec in enumerate(d_vec):
                f += self.morse_force(d_mag[count], *param[0:3]) * vec / d_mag[count] + \
                     self.bonding_force(d_mag[count], *param[7:11]) * vec / d_mag[count]
            dict_morse[atom_tag] = np.linalg.norm(f)

        for atom_tag, angle_list in dict_EF_angles.items():
            angle, f = [], 0
            if 'Cu' in atom_tag:
                my_param = param[3:5]
            else:
                my_param = param[5:7]
            for count, ang in enumerate(angle_list):
                f += self.harmonic(ang, *my_param)
            dict_EF_predicted_forces[atom_tag] = dict_morse[atom_tag] + f

            if 'O' in atom_tag:
                for count, ang in enumerate(angle_list):
                    f += self.harmonic(ang, *param[3:5])
            dict_EF_predicted_forces[atom_tag] = dict_morse[atom_tag] + f

        return dict_EF_predicted_forces
    """


if __name__ == '__main__':
    traj = read('/Users/jiaweiguo/Box/ECH289_Project/MFI.traj', '4:5')
    for atoms in traj:
        try:
            EF_analyzer = ExtraFrameworkAnalyzer(atoms)
            EF_analyzer.get_extraframework()
            print('atom index: \n', EF_analyzer.dict_EF)
            print('bonds: \n', EF_analyzer.get_all_bonds())
            print('angles: \n', EF_analyzer.get_all_angles())
            # print(EF_analyzer.get_angle_force_dir())
            print('DFT forces: \n', EF_analyzer.get_forces())
        except:
            print('Error')
