from maze.extra_framework_maker import ExtraFrameworkMaker
from maze.zeolite import PerfectZeolite, Zeolite
from ase.neighborlist import natural_cutoffs, NeighborList
from ase.io import write, read
from ase.visualize import view
from ase import Atoms
import os
from pathlib import Path


def main():
    sample_zeolite = 'AEI'
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

    output_dir = '/Users/jiaweiguo/Desktop'
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
            write(my_path, atoms, format='vasp')

    for site_name, traj in dict_2Al_ZH.items():
        output_dir2 = os.path.join(output_dir, sample_zeolite, '2Al', site_name)
        Path(output_dir2).mkdir(parents=True, exist_ok=True)
        for atoms in traj:
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = nl.get_neighbors(index_H[0])[0]
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val != index_H][0]
            output_filename = sample_zeolite + '_' + '2Al_' + site_name + '_O_' + str(neigh_o_indices[0]) + '_O_' + str(neigh_o_indices[1])
            my_path = os.path.join(output_dir2, output_filename)
            write(my_path, atoms, format='vasp')
    print('File writing is finished!')
    # ADD TAGS FOR O AS WELL (BESIDES INDEX)


if __name__ == '__main__':
    main()


