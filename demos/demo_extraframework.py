from maze.extra_framework_maker import ExtraFrameworkMaker
from ase.io import write, read
from ase.visualize import view
from ase import Atoms


def main():
    sample_zeolite = 'AEI'
    EFzeolite = ExtraFrameworkMaker(iza_code=sample_zeolite)
    # 1Al and 2Al replacement
    EFzeolite.make_extra_frameworks(print_statement=True)
    print(EFzeolite.t_site_indices)
    # print(len(EFzeolite.traj_1Al))
    print('Number of Unique Al pairs: ', len(EFzeolite.traj_2Al))
    # view(EFzeolite.traj_1Al)
    # view(EFzeolite.traj_2Al)

    # extra-framework insertion (extra-framework atoms input can be ASE Atoms object or user-defined traj files)
    zeolite = EFzeolite.traj_2Al[0]
    EF_atoms = Atoms('Co', positions=[[0, 0, 0]])
    # EF_atoms = read('/Users/jiaweiguo/Desktop/MAZE-sim-master/demos/CuOCu_cluster.traj', '0')
    inserted_atoms = EFzeolite.insert_ExtraFrameworkAtoms(zeolite, EF_atoms)
    view(inserted_atoms)
    print('Extra-framework atoms insertion is done!')

    # sample all Z-TM sites
    dict_Z_TM = EFzeolite.get_all_Z_TM(d_Z_TM=2.6, TM_type='Cu')  # Al-Cu distance
    all_traj = []
    for site_name, traj in dict_Z_TM.items():
        for each_traj in traj:
            all_traj.append(each_traj)
    view(all_traj)
    print('All Z-TM site sampling is done!')

    # sample all Bronsted sites
    traj_ZH = EFzeolite.get_all_Bronsted_sites()
    view(traj_ZH)
    print('All Bronsted site sampling is done!')


if __name__ == '__main__':
    main()
