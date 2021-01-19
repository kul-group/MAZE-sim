from ase.io import read, write
from ase import Atoms
from maze import Zeotype, OpenDefect, ImperfectZeotype
import os
from ase.visualize import view
from pathlib import Path
from collections import defaultdict
from typing import List
import numpy as np


# %%
def defect_maker(cif_dir, zeolite_code, output_dir, savefiles=True):
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    open_defect = OpenDefect(zeolite)
    unique_t_site_indices = {}
    for site_name, value in zeolite.site_to_atom_indices.items():
        if 'T' in site_name:
            unique_t_site_indices[site_name] = value[1]

    # build dictionary of open defects
    unique_t_site_to_od = defaultdict(list)
    for site_name, t_site in unique_t_site_indices.items():
        for o_index in open_defect.neighbor_list.get_neighbors(t_site)[0]:
            for si_index in open_defect.neighbor_list.get_neighbors(o_index)[0]:
                if si_index == t_site:
                    continue
                pos = open_defect.get_positions()[si_index]
                new_od = open_defect.delete_atoms([si_index]).cap_atoms()
                unique_t_site_to_od[site_name].append(new_od)

    # save T sites
    if savefiles:  #I made this change
        for site_name, od_list in unique_t_site_to_od.items():
            output_dir2 = os.path.join(output_dir, zeolite_code, site_name)
            Path(output_dir2).mkdir(parents=True, exist_ok=True)
            for index, od in enumerate(od_list):
                output_filename = zeolite_code + '_' + site_name + '_' + str(index) + '.traj'
                my_path = os.path.join(output_dir2, output_filename)
                write(my_path, od)

    return unique_t_site_to_od  # and most importantly I made this chnage


# helper functions for second defect function
def find_removed_atoms(iz: ImperfectZeotype) -> List[int]:
    """
    Finds the atoms removed from the iz by comparing it with parent MAZE-sim
    :param iz: imperfect MAZE-sim to check for missing atoms
    :return: list of indices of missing atoms (using parent indexing)
    """
    missing_list = []
    for atom in iz.parent_zeotype:
        value = iz.index_mapper.get_index(iz.parent_zeotype.name, iz.name, atom.index)
        if value is None:
            missing_list.append(atom.index)

    return missing_list

def find_neighbor_si(z: Zeotype, first_si_index: int):
    """
    Finds the first neighboring Si
    :param z: the MAZE-sim object
    :param first_si_index: the first Si index to find a neighbor too
    :return: the first neighboring si index
    """
    z.update_nl()
    for o_index in z.neighbor_list.get_neighbors(first_si_index)[0]:
        for si_index in z.neighbor_list.get_neighbors(o_index)[0]:
            if si_index == first_si_index:
                continue
            return si_index

def find_index_common_oxygen(iz, site_1: int, site_2: int) -> int:
    """
    Finds a common oxygen, if it exists, between two T sites
    :param iz: imperfect MAZE-sim (or subclass) containing T sites
    :param site_1: index of T site 1
    :param site_2: index of T site 2
    :return:
    """
    iz.update_nl()
    nl1 = iz.neighbor_list.get_neighbors(site_1)[0]
    nl2 = iz.neighbor_list.get_neighbors(site_2)[0]
    # find common oxygen
    for i in nl1:
        if i in nl2:
            if iz[i].symbol == 'O':
                return i

    assert False, 'No middle oxygen found!!'

def remove_two_T_sites(iz, site_1: int, site_2: int) -> ImperfectZeotype:
    """
    Removes two T sites that are adjacent to eachother
    :param iz: Impefect MAZE-sim with two T sites
    :param site_1: the index of the first site to remove
    :param site_2: the index of the second site to remove
    :return:
    """
    indices_to_remove = [site_1, site_2, find_index_common_oxygen(iz, site_1, site_2)]
    return iz.delete_atoms(indices_to_remove)

def second_defect(od_dict, T_site_str, list_pos, cif_name, out_path, savefile: bool =True):
    """
    Second defect creator
    :param od_dict: dictionary of open defects (this is to allow it to integrate with the other code better)
    :param T_site_str: The key for the open defects dictionary
    :param list_pos: The position in the list of od_dict[T_site_str]
    :param cif_name: the name of the cif file used to build the zeolite
    :param out_path: the output path
    :param savefile: bool if a file should be saved or not
    :return: an open defect with an additional adjacent site removed
    """
    od = od_dict[T_site_str][list_pos]  # get specific open defect
    complete_od = OpenDefect(od.parent_zeotype)  # create an opendefect object without the removed index
    removed_atom_index = find_removed_atoms(od)[0]

    neighbor_si_index = find_neighbor_si(complete_od, removed_atom_index)
    new_od = remove_two_T_sites(complete_od, removed_atom_index, neighbor_si_index).cap_atoms()
    my_path = os.path.join(out_path, cif_name + '_' + str(neighbor_si_index) + '.traj')
    if savefile:
        write(my_path, new_od)
    return new_od


# .add_atoms(Atoms('Ir', positions=[pos]), 'Ir')

# %%'

# %%


def insert_HM(cif_dir, out_path, cif_name, si_index):
    atoms = read(cif_dir + str(si_index) + '.traj', '-1')
    zeolite = Zeotype(atoms)
    open_defect = OpenDefect(zeolite)
    atoms_H = [a.index for a in open_defect if a.symbol in ['H']]
    pos_H = open_defect.get_positions()[atoms_H]
    cell_dim = atoms.get_cell_lengths_and_angles()[0:3]
    for row in pos_H:
        for index, item in enumerate(row):
            if item >= cell_dim[index]:
                item -= cell_dim[index]
                row[index] = item
            elif item < 0:
                item += cell_dim[index]
                row[index] = item
            else:
                continue
    pos = np.mean(pos_H, axis=0)
    new_od = open_defect.add_atoms(Atoms('Ir', positions=[pos]), 'Ir')
    # my_path = os.path.join(out_path, cif_name + '_' + str(si_index) + '_Ir.traj')
    # write(my_path, new_od)
    return np.mean(pos)


def insert_HM_2(open_defect, si_index):
    atoms_types, _ = open_defect.count_elements
    atoms_H = atoms_types['H']
    pos_H = open_defect.get_positions()[atoms_H]
    cell_dim = open_defect.get_cell_lengths_and_angles()[0:3]
    for row in pos_H:
        for index, item in enumerate(row):
            if item >= cell_dim[index]:
                item -= cell_dim[index]
                row[index] = item
            elif item < 0:
                item += cell_dim[index]
                row[index] = item
            else:
                continue
    pos = np.mean(pos_H, axis=0)
    new_od = open_defect.add_atoms(Atoms('Ir', positions=[pos]), 'Ir')
    return np.mean(pos)


if __name__ == "__main__":
    #defect_maker('/Users/jiaweiguo/Desktop/0125Proposal/BEA.cif', 'BEA', '/Users/jiaweiguo/Desktop/0125Proposal')
    od_dict = defect_maker('/data/BEA.cif', 'BEA', '//data/test_output',
                           savefiles=False)
    my_od = second_defect(od_dict, 'T3', 3, 'BEA', '//data/test_output', savefile=False)
    view(my_od)
    # second_defect(cif_dir, out_path, 'BEA_T1_3', 189)
    # second_defect(cif_dir, out_path, 'BEA_T1_3', 141)
    # second_defect(cif_dir, out_path, 'BEA_T1_3', 177)
    # cif_dir = '/Users/jiaweiguo/Desktop/0125Proposal/BEA/T1/BEA_T1_3.traj'
    # out_path = '/Users/jiaweiguo/Desktop/0125Proposal/BEA/T1/'
    # cif_dir = '/Users/jiaweiguo/Desktop/0125Proposal/BEA/T1/BEA_T1_3_'
    # out_path = '/Users/jiaweiguo/Desktop/0125Proposal/BEA/T1/'
    # insert_HM(cif_dir, out_path, 'BEA_T1_3', 141)
    # insert_HM(cif_dir, out_path, 'BEA_T1_3', 189)
    # insert_HM(cif_dir, out_path, 'BEA_T1_3', 177)
    #
