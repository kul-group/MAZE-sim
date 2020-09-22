from source.zeotype import Zeotype


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    atom_types = zeolite.get_atom_types()
    for key, val in atom_types.items():
        print(key, val)

    indices, count = zeolite.count_elements()
    print('count', dict(count))

if __name__ == '__main__':
    main()
