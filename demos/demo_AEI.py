from maze.extra_framework_maker import ExtraFrameworkMaker
from maze.zeolite import PerfectZeolite, Zeolite
from ase.neighborlist import natural_cutoffs, NeighborList
from ase.io import write, read
from ase.visualize import view
from ase import Atoms
import os
from pathlib import Path


def Bronsted_H_structures_reduced_AEI():
    sample_zeolite = 'AEI'
    AEI_path = '/Users/jiaweiguo/Desktop/AEI/AEI_reduced_random_T.cif'
    EFzeolite = ExtraFrameworkMaker('AEI', user_input_path=AEI_path)
    # print(EFzeolite.EFzeolite.site_to_atom_indices)

    for atom in EFzeolite.EFzeolite:
        if atom.symbol != 'O':
            atom.symbol = 'Si'

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

    output_dir0 = '/Users/jiaweiguo/Desktop/AEI'
    write(output_dir0 + '/reduced_AEI.vasp', EFzeolite.EFzeolite, format='vasp')

    output_dir1 = os.path.join(output_dir0, '1Al_H')
    Path(output_dir1).mkdir(parents=True, exist_ok=True)
    for site_name, traj in dict_1Al_ZH.items():
        for count, atoms in enumerate(traj):
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = nl.get_neighbors(index_H[0])[0]
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val != index_H][0]
            output_filename = sample_zeolite + '_' + '1Al_' + site_name + '_H_O' + str(neigh_o_indices)
            my_path = os.path.join(output_dir1, output_filename)
            write(my_path+'.vasp', atoms, format='vasp')

    output_dir2 = os.path.join(output_dir0, '2Al_H')
    Path(output_dir2).mkdir(parents=True, exist_ok=True)
    for site_name, traj in dict_2Al_ZH.items():
        for atoms in traj:
            index_H = [a.index for a in atoms if a.symbol == 'H']
            index_O = [a.index for a in atoms if a.symbol == 'O']
            nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
            nl.update(atoms)
            neigh_indices = list(nl.get_neighbors(index_H[0])[0]) + list(nl.get_neighbors(index_H[1])[0])
            neigh_o_indices = [val for val in neigh_indices if val in index_O and val not in index_H]
            assert len(neigh_o_indices) == 2
            output_filename = sample_zeolite + '_' + '2Al_' + site_name + '_O_' + str(neigh_o_indices[0]) + '_O_' + str(neigh_o_indices[1])
            my_path = os.path.join(output_dir2, output_filename)
            write(my_path+'.vasp', atoms, format='vasp')

    print('File writing is finished!')


def OSDA_structures_reduced_AEI():
    sample_zeolite = 'AEI'
    file = 'cis+transCuSSZ39_POSCAR'
    atoms = read('/Users/jiaweiguo/Desktop/AEI/%s.vasp' % file)

    """
    AEI_path = '/Users/jiaweiguo/Desktop/AEI/cis+transCuSSZ39_POSCAR.vasp'

    EFzeolite = ExtraFrameworkMaker('AEI', user_input_path=AEI_path)
    # print(EFzeolite.EFzeolite.site_to_atom_indices)

    for atom in EFzeolite.EFzeolite:
        if atom.symbol != 'O':
            atom.symbol = 'Si'

    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))
    """


if __name__ == '__main__':
    # main()
    sample_zeolite = 'AEI'
    file = 'CU-SSZ39_trans_with_Na' # 'CU-SSZ39_cis_with_Na'
    # 'doubletrans_Cu-SSZ39_POSCAR' # 'doublecis_CuSSZ39_POSCAR' # 'cis+transCuSSZ39_POSCAR'

    atoms = read('/Users/jiaweiguo/Desktop/AEI/%s.cif' % file) # .vasp
    view(atoms)

    cluster_index = [atom.index for atom in atoms if atom.symbol not in ['Si', 'O', 'Al']]
    cluster = atoms[cluster_index]
    # view(cluster)

    AEI_path = '/Users/jiaweiguo/Desktop/AEI/AEI_reduced_random_T.cif'
    EFzeolite = ExtraFrameworkMaker('AEI', user_input_path=AEI_path)
    # print(EFzeolite.EFzeolite.site_to_atom_indices)

    for atom in EFzeolite.EFzeolite:
        if atom.symbol != 'O':
            atom.symbol = 'Si'

    EFzeolite.make_extra_frameworks(replace_1Al=True, replace_2Al=True, print_statement=True)
    print(EFzeolite.t_site_indices)
    print(EFzeolite.t_site_indices_count)
    print('Number of unique Al pairs: ', len(EFzeolite.traj_2Al))

    output_dir0 = '/Users/jiaweiguo/Desktop/AEI'
    output_dir1 = os.path.join(output_dir0, '2Al_%s' % file)
    Path(output_dir1).mkdir(parents=True, exist_ok=True)
    for site_name, atoms in EFzeolite.dict_2Al_replaced.items():
        atoms = atoms + cluster
        output_filename = sample_zeolite + '_' + '2Al_' + site_name
        my_path = os.path.join(output_dir1, output_filename)
        write(my_path + '.vasp', atoms, format='vasp')
