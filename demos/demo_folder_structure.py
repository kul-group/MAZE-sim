# generating folder structures for all zeolites

from maze.extra_framework_maker import ExtraFrameworkMaker
from ase.neighborlist import natural_cutoffs, NeighborList
import os
from pathlib import Path
from ase.io import write, read
from ase.visualize import view


def save_zeolites(sample_zeolite: str):
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))

    # sample all Bronsted sites
    dict_1Al_ZH = EFzeolite.get_all_Bronsted_sites(case_1Al=True)
    print('Single Al Bronsted site sampling is done!')
    dict_2Al_ZH = EFzeolite.get_all_Bronsted_sites(case_2Al=True)
    print('Two Al Bronsted site sampling is done!')

    output_dir = '/Users/jiaweiguo/Box/all_zeo_database'

    for site_name, traj in dict_1Al_ZH.items():
        output_dir1 = os.path.join(output_dir, sample_zeolite, '1Al', site_name, 'H')
        Path(output_dir1).mkdir(parents=True, exist_ok=True)

        for count, atoms in enumerate(traj):
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = nl.get_neighbors(index_H[0])[0]
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val != index_H][0]
            output_filename = sample_zeolite + '_' + '1Al_' + site_name + '_H_O' + str(neigh_o_indices)
            # use O tag instead of indices (maybe)
            my_path = os.path.join(output_dir1, output_filename)
            write(my_path+'.traj', atoms)

    for site_name, traj in dict_2Al_ZH.items():
        output_dir2 = os.path.join(output_dir, sample_zeolite, '2Al', site_name, 'H')
        Path(output_dir2).mkdir(parents=True, exist_ok=True)

        for atoms in traj:
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = list(nl.get_neighbors(index_H[0])[0]) + list(nl.get_neighbors(index_H[1])[0])
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val not in index_H]
            assert len(neigh_o_indices) == 2
            output_filename = sample_zeolite + '_' + '2Al_H_' + site_name + '_O' + str(neigh_o_indices[0]) + '_O' + str(neigh_o_indices[1])
            my_path = os.path.join(output_dir2, output_filename)
            # write(my_path, atoms, format='vasp') write as poscar
            write(my_path+'.traj', atoms)

    # insert extra-framework CuOCu
    EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
    EF_type = 'CuOCu'
    for site_name, zeolite in EFzeolite.dict_2Al_replaced.items():
        output_dir3 = os.path.join(output_dir, sample_zeolite, '2Al', site_name, EF_type)
        Path(output_dir3).mkdir(parents=True, exist_ok=True)

        atoms = EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms)
        output_filename = sample_zeolite + '_' + '2Al_CuOCu_' + site_name
        my_path = os.path.join(output_dir3, output_filename)
        write(my_path+'.traj', atoms)

    print('Extra-framework atoms insertion is done!')

    print('File writing is finished!')
    # ase.build('CuOCu')
    #TODO: ADD TAGS FOR O AS WELL (BESIDES INDEX)
    #TODO: ADD FOLDER STRUCTURE FUNCTION IN


def save_zeolites_v2(sample_zeolite: str, output_dir):

    # save bare zeolite structure
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)
    output_dir0 = os.path.join(output_dir, sample_zeolite)
    Path(output_dir0).mkdir(parents=True, exist_ok=True)
    my_path = os.path.join(output_dir0, sample_zeolite)
    write(my_path + '.traj', EFzeolite.EFzeolite)

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))

    for site_name, traj in EFzeolite.dict_1Al_replaced.items():
        output_dir1 = os.path.join(output_dir, sample_zeolite, '1Al', site_name)
        Path(output_dir1).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir1, sample_zeolite + '_' + '1Al_' + site_name)
        write(my_path + '.traj', EFzeolite.EFzeolite)

    for site_name, traj in EFzeolite.dict_2Al_replaced.items():
        output_dir2 = os.path.join(output_dir, sample_zeolite, '2Al', site_name)
        Path(output_dir2).mkdir(parents=True, exist_ok=True)
        my_path = os.path.join(output_dir2, sample_zeolite + '_' + '1Al_' + site_name)
        write(my_path + '.traj', EFzeolite.EFzeolite)

    """
    # sample all Bronsted sites
    dict_1Al_ZH = EFzeolite.get_all_Bronsted_sites(case_1Al=True)
    print('Single Al Bronsted site sampling is done!')
    dict_2Al_ZH = EFzeolite.get_all_Bronsted_sites(case_2Al=True)
    print('Two Al Bronsted site sampling is done!')


    for site_name, traj in dict_1Al_ZH.items():
        output_dir1 = os.path.join(output_dir, sample_zeolite, '1Al', site_name, 'H')
        Path(output_dir1).mkdir(parents=True, exist_ok=True)

        for count, atoms in enumerate(traj):
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = nl.get_neighbors(index_H[0])[0]
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val != index_H][0]
            output_filename = sample_zeolite + '_' + '1Al_' + site_name + '_H_O' + str(neigh_o_indices)
            # use O tag instead of indices (maybe)
            my_path = os.path.join(output_dir1, output_filename)
            write(my_path+'.traj', atoms)

    for site_name, traj in dict_2Al_ZH.items():
        output_dir2 = os.path.join(output_dir, sample_zeolite, '2Al', site_name, 'H')
        Path(output_dir2).mkdir(parents=True, exist_ok=True)

        for atoms in traj:
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = list(nl.get_neighbors(index_H[0])[0]) + list(nl.get_neighbors(index_H[1])[0])
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val not in index_H]
            assert len(neigh_o_indices) == 2
            output_filename = sample_zeolite + '_' + '2Al_H_' + site_name + '_O' + str(neigh_o_indices[0]) + '_O' + str(neigh_o_indices[1])
            my_path = os.path.join(output_dir2, output_filename)
            # write(my_path, atoms, format='vasp') write as poscar
            write(my_path+'.traj', atoms)

    # insert extra-framework CuOCu
    EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
    EF_type = 'CuOCu'
    for site_name, zeolite in EFzeolite.dict_2Al_replaced.items():
        output_dir3 = os.path.join(output_dir, sample_zeolite, '2Al', site_name, EF_type)
        Path(output_dir3).mkdir(parents=True, exist_ok=True)

        atoms = EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms)
        output_filename = sample_zeolite + '_' + '2Al_CuOCu_' + site_name
        my_path = os.path.join(output_dir3, output_filename)
        write(my_path+'.traj', atoms)

    print('Extra-framework atoms insertion is done!')

    print('File writing is finished!')
    """


if __name__ == '__main__':
    # save_zeolites('CHA')
    save_zeolites_v2('CHA', '/Users/jiaweiguo/Box/all_zeo_database')
