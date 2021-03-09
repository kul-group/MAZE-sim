from maze.extra_framework import ExtraFramework
from ase.io import write, read
from ase.visualize import view


def main():
    cif_dir = 'MFI.cif'
    my_zeolite = ExtraFramework(cif_dir=cif_dir)
    # unique_t_site_indices = my_zeolite.get_unique_t_sites(cif_dir)
    my_zeolite.create_1Al_replacement()
    my_zeolite.create_2Al_replacement()
    # write('MFI_2Al_replaced.traj', my_zeolite.traj_2Al)
    print(my_zeolite.count)

    traj = read('MFI_2Al_replaced.traj', ':')
    traj_CuOCu = []
    for atoms in traj:
        traj_CuOCu.append(my_zeolite.insert_CuOCu(atoms))
    view(traj_CuOCu)


if __name__ == '__main__':
    main()
