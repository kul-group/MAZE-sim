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
from statistics import mode


class ExtraFrameworkAnalyzer(object):

    def __init__(self, atoms):
        self.TM_list = ['Pt', 'Cu', 'Co', 'Pd', 'Fe', 'Cr', 'Rh', 'Ru']
        self.dict_EF = {}
        self.atoms = PerfectZeolite(atoms)

    def get_extraframework_cluster(self, predefined_centering_o=None):
        """
        extract extra-framework cluster including Cu-O-Cu, 2 Al, and 8 O around the Als (13 atoms total)
        :param predefined_centering_o: get the mode of all possible centering O index by training a bunch of
        configurations for the same zeolite, to remove bias.
        """
        index_EF_TM = [a.index for a in self.atoms if a.symbol in self.TM_list]
        index_Al = [a.index for a in self.atoms if a.symbol == 'Al']

        Al_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_Al[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_Al[1])[0]))
        Al_neigh_list = [x for x in Al_neigh_list if self.atoms[x].symbol == 'O']

        if predefined_centering_o is not None:
            centering_o = copy.copy(predefined_centering_o)
        else:
            TM_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_EF_TM[0])[0],
                                            self.atoms.neighbor_list.get_neighbors(index_EF_TM[1])[0]))
            centering_o = [[x for x in TM_neigh_list if list(TM_neigh_list).count(x) > 1 and x not in Al_neigh_list][0]]

        assert len(index_EF_TM) == 2
        assert len(index_Al) == 2
        assert len(centering_o) == 1
        assert len(Al_neigh_list) == 8

        return Al_neigh_list + index_Al + index_EF_TM + centering_o


def check_centering_o_index():
    zeo_list = ['MWW', 'CHA', 'AEI', 'BEA', 'LTA', 'MAZ', 'MFI', 'MOR', 'RHO', 'SOD']

    for zeolite in zeo_list:
        traj = read('/Users/jiaweiguo/Box/openMM_FF/%s_minE.traj' % zeolite, ':')
        # view(traj)
        EF_O_index, count = [], 0
        for atoms in traj:
            try:
                EFAnalyzer = ExtraFrameworkAnalyzer(atoms)
                EF_O_index.append(EFAnalyzer.get_extraframework_cluster()[-1])
                count += 1
            except:
                ...
        unique_index = np.unique(EF_O_index)
        print(zeolite, count, list(unique_index),
              [EF_O_index.count(unique_index[count]) for count in range(len(unique_index))])
        # index check, if len(np.unique(EF_O_index)) > 1, remove


if __name__ == '__main__':
    atoms = read('/Users/jiaweiguo/Box/openMM_FF/%s_minE.traj' %'CHA', '0')
    EFAnalyzer = ExtraFrameworkAnalyzer(atoms)
    print(EFAnalyzer.get_extraframework_cluster([144])[-1])
    
