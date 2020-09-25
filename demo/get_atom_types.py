from source.zeotype import Zeotype


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    print('t_site_to_atom_indices', zeolite.t_site_to_atom_indices)
    print('atom_indices_to_t_site', zeolite.atom_indices_to_t_site)
    atom_types = zeolite.get_atom_types()
    for key, val in atom_types.items():
        print(key, val)

    indices, count = zeolite.count_elements()
    print('atom type count', dict(count))

if __name__ == '__main__':
    main()
