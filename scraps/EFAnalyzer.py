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

    def get_extraframework_atoms(self):
        """
        extract extra-framework cluster including Cu-O-Cu, 2 Al, and 8 O around the Als (13 atoms total)
        """
        index_EF_TM = [a.index for a in self.atoms if a.symbol in self.TM_list]
        index_Al = [a.index for a in self.atoms if a.symbol == 'Al']
        
        TM_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_EF_TM[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_EF_TM[1])[0]))

        # This might be wrong! Not always return desired O index!
        # works for All CHA structures at least
        centering_o = [[x for x in TM_neigh_list if list(TM_neigh_list).count(x) > 1][0]]
        # print(centering_o)

        self.atoms.update_nl(1.2)
        Al_neigh_list = np.concatenate((self.atoms.neighbor_list.get_neighbors(index_Al[0])[0],
                                        self.atoms.neighbor_list.get_neighbors(index_Al[1])[0]))
        Al_neigh_list = [x for x in Al_neigh_list if self.atoms[x].symbol == 'O']
        
        assert len(index_EF_TM) == 2
        assert len(index_Al) == 2
        assert len(centering_o) == 1
        assert len(Al_neigh_list) == 8

        return Al_neigh_list + index_Al + index_EF_TM + centering_o


if __name__ == '__main__':
    traj = read('/Users/jiaweiguo/Box/openMM_FF/CHA_minE.traj', ':')
    for atoms in traj:
        EFAnalyzer = ExtraFrameworkAnalyzer(atoms)
        print(EFAnalyzer.get_extraframework_atoms())
    # view(traj)
