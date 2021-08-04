from maze.zeolite import PerfectZeolite, Zeolite
from maze.io_zeolite import read_vasp
from ase import Atoms, db
from typing import Union, Tuple
from collections import defaultdict
from ase.neighborlist import natural_cutoffs, NeighborList, mic
from ase import Atoms
import numpy as np
from ase.visualize import view
import copy as copy
from ase.io import write, read
import itertools


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

    @staticmethod
    def recentering_atoms(atoms, pos_to_center):
        """ This function recenters the atoms object by translating the the input position "pos_to_center" to the center
        of the cell.
        :param atoms: zeolite backbone with two Al inserted
        :param pos_to_center: some positions to be centered
        :return: recentered atoms object, translation vector
        """
        vec_translate = np.matmul([0.5, 0.5, 0.5], atoms.get_cell()) - pos_to_center
        atoms.translate(vec_translate)
        atoms.wrap()
        return atoms, vec_translate

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
    """
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
        centering_o = [288] #[[x for x in TM_neigh_list if list(TM_neigh_list).count(x) > 1][0]]
        self.centering_atom_index = centering_o[0]
        assert len(centering_o) == 1
        assert self.atoms[centering_o].symbols == 'O'

        self.EF_indices = list(centering_o)
        self.EF_indices.extend([value for value in index_EF_TM])

        return np.unique(list(Al_neigh_list) + centering_o + index_Al + index_EF_TM)
